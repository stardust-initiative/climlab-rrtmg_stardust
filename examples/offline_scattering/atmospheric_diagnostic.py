"""Export clear-sky native operators and compare with independent offline algebra."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

# Set before importing the backend. Dump writes are deliberately serial.
if os.environ.get('OMP_NUM_THREADS', '1') != '1':
    raise RuntimeError('Operator export requires OMP_NUM_THREADS=1')
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np

from offline_scattering import Slab, reflected, solve_interfaces, stack

FORMAT = 17001


def read_dump(path):
    """Read versioned little-endian gfortran stream; reject incomplete records."""
    records = []
    with Path(path).open('rb') as stream:
        def read(dtype, count):
            result = np.fromfile(stream, dtype=dtype, count=count)
            if result.size != count:
                raise ValueError('Truncated operator record')
            return result
        while True:
            header = np.fromfile(stream, dtype='<i4', count=4)
            if not len(header):
                break
            if len(header) != 4 or header[0] != FORMAT:
                raise ValueError('Unsupported operator format')
            _, n, mode, ng = map(int, header)
            if not (1 <= n <= 1000 and ng == 112 and mode in (0, 1)):
                raise ValueError('Unexpected dimensions or mode')
            mu = float(read('<f8', 1)[0])
            data, aero, bands, weights = [], [], [], []
            for g in range(ng):
                iw, band = read('<i4', 2)
                if iw != g+1 or not 1 <= band <= 14:
                    raise ValueError('Invalid spectral ordering')
                weights.append(float(read('<f8', 1)[0]))
                data.append(read('<f8', 10*(n+1)).reshape(10, n+1))
                aero.append(read('<f8', 4*n).reshape(4, n))
                bands.append(int(band))
            d, a, w = np.array(data), np.array(aero), np.array(weights)
            if not (np.isfinite(d).all() and np.isfinite(a).all() and np.isfinite(w).all()
                    and np.all(w >= 0) and 0 < mu <= 1):
                raise ValueError('Invalid operator values')
            records.append(dict(mode=mode, mu=mu, data=d, aerosol=a, bands=np.array(bands), weights=w))
    if not records:
        raise ValueError('Empty operator dump')
    return records


def evaluate(record, first, stop):
    """Derive lower Ab/Ad and full fluxes using exactly the exported beam."""
    d, w = record['data'], record['weights']
    ng, _, interfaces = d.shape
    n = interfaces-1
    if not 0 <= first < stop <= n:
        raise ValueError('Invalid aerosol-region boundary')
    ab, ad, upward, downward = [], [], [], []
    interface_error = top_error = lower_error = 0.
    for g in range(ng):
        cells = [Slab.symmetric_cell(*d[g, :5, k]) for k in range(n)]
        # Last native coefficient row is the lower surface response.
        surface = Slab(d[g, 0, n], d[g, 1, n], 0., 0., 0., 0., 0.)
        lower = stack(cells[stop:] + [surface])
        ab.append(lower.reflect_beam); ad.append(lower.reflect_down)
        lower_error = max(lower_error, abs(ab[-1]-d[g, 6, stop]), abs(ad[-1]-d[g, 7, stop]))
        combined = stack([stack(cells[:first]), stack(cells[first:stop])])
        top = reflected(combined, ab[-1], ad[-1])
        top_error = max(top_error, abs(top-d[g, 9, 0]))
        beam, diffuse, up = solve_interfaces(cells + [surface])
        down = (beam+diffuse)[:-1]; up = up[:-1]
        interface_error = max(interface_error, float(np.max(abs(up-d[g, 9]))), float(np.max(abs(down-d[g, 8]))))
        upward.append(up); downward.append(down)
    up, down = np.array(upward), np.array(downward)
    if max(interface_error, top_error, lower_error) > 2e-11:
        raise AssertionError('Offline/native coefficient comparison failed')
    return dict(lower_direct=np.array(ab), lower_diffuse=np.array(ad), offline_up=up, offline_down=down,
                native_up=d[:, 9], native_down=d[:, 8], weights=w, bands=record['bands'],
                max_interface_error=interface_error, max_top_error=top_error, max_lower_response_error=lower_error,
                broadband_up=w@up, broadband_down=w@down)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='new output directory')
    parser.add_argument('--build-manifest', type=Path, required=True, help='diagnostic_build.json from source preparation')
    args = parser.parse_args()
    manifest = json.loads(args.build_manifest.read_text())
    if manifest.get('format_version') != FORMAT:
        parser.error('Build manifest format differs from exporter')
    if sys.byteorder != 'little':
        parser.error('The documented diagnostic stream currently requires little-endian GNU Fortran')
    args.output.mkdir(parents=True, exist_ok=False)
    from atmosphere_fixture import column, pressure_bounds
    from climlab_rrtmg import rrtmg_sw
    import climlab_rrtmg.rrtmg_sw._rrtmg_sw as native
    outputs, records, raw_records = {}, {}, {}
    layers = [(0., 1., 0., 1.)]*6
    # Effective coefficients are already top-to-bottom; atmosphere arrays reverse.
    layers[3] = (.08, .88, .12, .84)
    layers[4] = (.04, .93, .07, .90)
    occupied = np.flatnonzero(np.any(np.array(layers) != [0., 1., 0., 1.], axis=1))
    first, stop = int(occupied[0]), int(occupied[-1]+1)
    for label, coefficients, enabled in [('clean', layers, False), ('aerosol', layers, True),
                                          ('identity', [(0., 1., 0., 1.)]*6, True)]:
        dump = args.output / (label+'.bin')
        os.environ['RRTMG_EXAMPLE_OPERATOR_DUMP'] = str(dump.resolve())
        native_result = column(coefficients, enabled=enabled, pressure=1000., mu=.5)
        if not all(np.isfinite(value).all() for value in native_result):
            raise AssertionError('Nonfinite native output, including heating; check initializer ABI')
        if not dump.exists():
            raise RuntimeError('Backend is not the separately instrumented diagnostic build')
        group = read_dump(dump)
        record = group[0]
        for other in group[1:]:
            if other['mode'] != record['mode'] or other['mu'] != record['mu']:
                raise AssertionError('Mixed native calls in one export')
            for key in ('data', 'aerosol', 'weights', 'bands'):
                np.testing.assert_array_equal(record[key], other[key])
        expected_mode = int(enabled)
        if record['mode'] != expected_mode:
            raise AssertionError('Unexpected native mode')
        if label == 'aerosol':
            expected = np.array(coefficients).T
            np.testing.assert_allclose(record['aerosol'], np.broadcast_to(expected, record['aerosol'].shape), atol=1e-14, rtol=0)
            active = np.flatnonzero(np.any(record['aerosol'] != np.array([0.,1.,0.,1.])[None,:,None], axis=(0,1)))
            first, stop = int(active[0]), int(active[-1]+1)
        result = evaluate(record, first, stop)
        # Restore physical weights and verify against full public wrapper outputs.
        np.testing.assert_allclose(result['broadband_up'], native_result[3][0, ::-1], atol=2e-9, rtol=1e-12)
        np.testing.assert_allclose(result['broadband_down'], native_result[4][0, ::-1], atol=2e-9, rtol=1e-12)
        np.savez(args.output/(label+'.npz'), **result, operators=record['data'], aerosol=record['aerosol'])
        outputs[label] = {k: v for k, v in result.items() if k.startswith('max_')}
        outputs[label]['dump_sha256'] = hashlib.sha256(dump.read_bytes()).hexdigest()
        outputs[label]['native_records'] = len(group)
        records[label] = result
        raw_records[label] = record
    np.testing.assert_allclose(records['clean']['broadband_up'], records['identity']['broadband_up'], atol=2e-10, rtol=0)
    np.testing.assert_allclose(records['clean']['broadband_down'], records['identity']['broadband_down'], atol=2e-10, rtol=0)
    np.testing.assert_allclose(raw_records['aerosol']['data'][:,4,:-1],
                               raw_records['clean']['data'][:,4,:-1] * raw_records['aerosol']['aerosol'][:,1,:],
                               atol=1e-14, rtol=1e-13)
    for key in ('lower_direct', 'lower_diffuse'):
        np.testing.assert_allclose(records['clean'][key], records['aerosol'][key], atol=1e-14, rtol=0)
    summary = dict(format_version=FORMAT, build_manifest=manifest,
                   build_manifest_sha256=hashlib.sha256(args.build_manifest.read_bytes()).hexdigest(), scope='clear-sky algebra with native gas/Rayleigh coefficients; not physical validation',
                   backend_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                   fixture_sha256=hashlib.sha256((Path(__file__).parent/'atmosphere_fixture.py').read_bytes()).hexdigest(),
                   pressure_bounds_hpa_top_to_bottom=pressure_bounds(len(layers))[::-1].tolist(),
                   aerosol_first_cell=first, aerosol_stop_cell_exclusive=stop, zenith_cosine=.5,
                   layers_order='top-to-bottom', spectral_weights_units='W m-2', coefficient_units='dimensionless',
                   zero_identity_max_error=float(np.max(abs(records['clean']['broadband_up']-records['identity']['broadband_up']))),
                   aerosol_toa_net_down_forcing_wm2=float(records['clean']['broadband_up'][0]-records['aerosol']['broadband_up'][0]),
                   results=outputs)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
