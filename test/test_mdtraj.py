import os
import unittest

import mdtraj


class MDTrajSugarTest(unittest.TestCase):
    gro_filename = 'test/data/sugar.gro'
    xtc_filename = 'test/data/sugar.xtc'
    n_atoms = 17

    topology = mdtraj.load(gro_filename).topology

    def test_estimate_frames(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            estimated_frames = xtc._estimate_n_frames_from_filesize(
                os.stat(self.xtc_filename).st_size)

            self.assertEqual(840, int(estimated_frames))

    def test_load_1_frame_manual(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            xyz, time, step, box, status = xtc._read(1, None, 1)

            self.assertEqual((1, self.n_atoms, 3), xyz.shape)
            self.assertEqual((1, 3, 3), box.shape)
            self.assertEqual((1,), time.shape)

    def test_load_1_frame(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            xyz, time, step, box = xtc.read(1)

            self.assertEqual((1, self.n_atoms, 3), xyz.shape)
            self.assertEqual((1, 3, 3), box.shape)
            self.assertEqual((1,), time.shape)

    def test_load_1000_frames(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            xyz, time, step, box = xtc.read(1000)

            self.assertEqual((1000, self.n_atoms, 3), xyz.shape)
            self.assertEqual((1000, 3, 3), box.shape)
            self.assertEqual((1000,), time.shape)

    def test_load_2000_frames(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            xyz, time, step, box = xtc.read(2000)

            self.assertEqual((1001, self.n_atoms, 3), xyz.shape)
            self.assertEqual((1001, 3, 3), box.shape)
            self.assertEqual((1001,), time.shape)

    def test_load_2000_frames_individually(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            for i in range(2000):
                xyz, time, step, box, status = xtc._read(1, None, 1)

                print(i)
                print(xyz)
                print(time)
                print(step)
                print(box)
                print(status)
                print('\n\n')

                # On failing macos build this succeeds for all 2000 frames
                # There are only 1000 frames in this trajctory...
                # These arrays should all have zero as the first dimension once we run out
                self.assertEqual((1, self.n_atoms, 3), xyz.shape)
                self.assertEqual((1, 3, 3), box.shape)
                self.assertEqual((1,), time.shape)

        self.fail()

    def test_load_xtc(self):
        with mdtraj.formats.XTCTrajectoryFile(self.xtc_filename) as xtc:
            traj = xtc.read_as_traj(self.topology)

            self.assertEqual(1001, traj.n_frames)


if __name__ == '__main__':
    unittest.main()
