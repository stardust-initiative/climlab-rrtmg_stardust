"""Prepare an opt-in diagnostic source copy; never edit the installed backend."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

DECLARATIONS = """      logical :: dump_active
      integer :: dump_unit, dump_status
      character(len=1024) :: dump_path

"""
INITIALIZE = """      call get_environment_variable('RRTMG_EXAMPLE_OPERATOR_DUMP',dump_path,status=dump_status)
      dump_active = dump_status == 0 .and. len_trim(dump_path) > 0
      if (dump_active) then
         if (any(pcldfmc > 0._rb)) error stop 'Diagnostic supports clear sky only'
         open(newunit=dump_unit,file=trim(dump_path),access='stream', &
              form='unformatted',status='unknown',position='append')
         write(dump_unit) 17001, nlayers, add_aero_layer, ngptsw, prmu0
      endif
"""
DUMP = """      if (dump_active) then
         write(dump_unit) iw, ibm, zincflx(iw)
         write(dump_unit) zrefc, zrefdc, ztrac, ztradc, zdbtc, ztdbtc
         write(dump_unit) zrupc, zrupdc, zcd(:,iw), zcu(:,iw)
         if (add_aero_layer == 1) then
            write(dump_unit) r_mu_col_g, t_mu_col_g, r_bar_col_g, t_bar_col_g
         else
            write(dump_unit) [(0._rb,jk=1,klev)], [(1._rb,jk=1,klev)]
            write(dump_unit) [(0._rb,jk=1,klev)], [(1._rb,jk=1,klev)]
         endif
      endif

"""

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="new directory, must not exist")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.destination = args.destination.resolve()
    if args.destination == root or root in args.destination.parents:
        parser.error("destination must be outside the source repository")
    if not args.destination.parent.is_dir():
        parser.error("destination parent must already exist")
    if args.destination.exists():
        parser.error("destination already exists; use a fresh directory")
    relative = Path("climlab_rrtmg/rrtmg_sw/sourcemods/rrtmg_sw_spcvmc.f90")
    text = (root / relative).read_text()
    for required in ["zdbtmc = zdbtmc * t_mu_col_g(jk)", "zdbtmo = zdbtmo * t_mu_col_g(jk)"]:
        if required not in text:
            raise RuntimeError("This diagnostic requires the corrected direct-beam source")
    original = hashlib.sha256(text.encode()).hexdigest()
    edits = [("      integer(kind=im)  :: klev", DECLARATIONS + "      integer(kind=im)  :: klev"),
             ("      ib1 = istart", INITIALIZE + "      ib1 = istart"),
             ("! Vertical quadrature for cloudy fluxes", DUMP + "! Vertical quadrature for cloudy fluxes"),
             ("      end subroutine spcvmc_sw", "      if (dump_active) close(dump_unit)\n      end subroutine spcvmc_sw")]
    for needle, replacement in edits:
        if text.count(needle) != 1:
            raise RuntimeError("Source layout changed; review diagnostic insertion: " + needle)
        text = text.replace(needle, replacement)
    shutil.copytree(root, args.destination, ignore=shutil.ignore_patterns(
        ".git", "build*", ".mesonpy-*", "__pycache__", "*.so", "*.pyc", ".pytest_cache"))
    source = args.destination / relative
    source.write_text(text)
    (args.destination / "diagnostic_build.json").write_text(json.dumps({
        "format_version": 17001, "scope": "clear-sky coefficients only; OMP_NUM_THREADS=1",
        "original_spcvmc_sha256": original,
        "driver_sha256": hashlib.sha256((root / "climlab_rrtmg/rrtmg_sw/Driver.f90").read_bytes()).hexdigest(),
        "diagnostic_spcvmc_sha256": hashlib.sha256(text.encode()).hexdigest()
    }, indent=2) + "\n")
    print(args.destination)

if __name__ == "__main__":
    main()
