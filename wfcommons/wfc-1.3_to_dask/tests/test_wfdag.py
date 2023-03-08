import unittest

from wfc2dask.wfctask import WFCTask
from wfc2dask.wfdag import WFDAG


import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


class TestWFDAG(unittest.TestCase):
    def test_sequence(self):
        in_fn = "samples/unittests/hello-world-sequence.json"
        tasks, wfname = WFCTask.load(in_fn)
        wfdag = WFDAG(wfname)
        for task in tasks:
            wfdag.add_task(task)
        wfdag.dask_codelines()

    def test_join(self):
        in_fn = "samples/unittests/hello-world-join.json"
        tasks, wfname = WFCTask.load(in_fn)
        wfdag = WFDAG(wfname)
        for task in tasks:
            wfdag.add_task(task)
        wfdag.dask_codelines()
        pass

    def test_big(self):
        in_fn = "samples/others/makeflow-instances/blast-chameleon-large-004.json"
        tasks, wfname = WFCTask.load(in_fn)
        wfdag = WFDAG(wfname)
        for task in tasks:
            wfdag.add_task(task)
        wfdag.dask_codelines()
        pass


if __name__ == '__main__':
    unittest.main()
