"""Small stream-format checks; native zero-identity is exercised by the workflow."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from atmospheric_diagnostic import FORMAT, read_dump


class StreamChecks(unittest.TestCase):
    def read(self, data):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'record.bin'
            path.write_bytes(data)
            return read_dump(path)

    def test_empty_rejected(self):
        with self.assertRaises(ValueError):
            self.read(b'')

    def test_unknown_format_rejected(self):
        with self.assertRaises(ValueError):
            self.read(np.array([1, 1, 0, 112], '<i4').tobytes())

    def test_truncation_rejected(self):
        with self.assertRaises(ValueError):
            self.read(np.array([FORMAT, 1, 0, 112], '<i4').tobytes())

    def test_valid_record(self):
        data = np.array([FORMAT, 1, 0, 112], '<i4').tobytes() + np.array([.5], '<f8').tobytes()
        for i in range(112):
            data += np.array([i+1, 1], '<i4').tobytes()
            data += np.array([1.], '<f8').tobytes()
            data += np.zeros(20, '<f8').tobytes() + np.zeros(4, '<f8').tobytes()
        records = self.read(data)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['data'].shape, (112, 10, 2))
        self.assertEqual(records[0]['weights'].sum(), 112.)


if __name__ == '__main__':
    unittest.main()
