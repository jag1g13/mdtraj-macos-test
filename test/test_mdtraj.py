import os
import unittest

import mdtraj


class MDTrajSugarTest(unittest.TestCase):
    gro_filename = 'test/data/sugar.gro'
    xtc_filename = 'test/data/sugar.xtc'

    topology = mdtraj.load(gro_filename).topology

    def test_estimate_frames(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as sugar_xtc:
            estimated_frames = sugar_xtc._estimate_n_frames_from_filesize(
                os.stat(self.xtc_filename).st_size)

            self.assertEqual(840, int(estimated_frames))

    def test_load_sugar_xtc(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as sugar_xtc:
            traj = sugar_xtc.read_as_traj(self.topology)

            self.assertEqual(1001, traj.n_frames)


if __name__ == '__main__':
    unittest.main()
