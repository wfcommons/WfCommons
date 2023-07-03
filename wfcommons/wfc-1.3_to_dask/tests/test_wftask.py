import unittest

from wfc2dask.wfctask import WFCTask


class TestWFTask(unittest.TestCase):
    def test_sequence(self):
        in_fn = "samples/unittests/hello-world-sequence.json"
        tasks = WFCTask.load(in_fn)
        pass

    def test_join(self):
        in_fn = "samples/unittests/hello-world-join.json"
        tasks = WFCTask.load(in_fn)
        pass


if __name__ == '__main__':
    unittest.main()
