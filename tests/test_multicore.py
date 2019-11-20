"""Test libs for multicore operations."""
import pytest
import multiprocessing
import subprocess

from idpconfgen.libs import libmulticore as LM

from . import tcommons

def test_Task_1():
    LM.Task()


class TestSubprocessTask:
    
    @pytest.mark.parametrize(
        'in1,in2,expected',
        [
            ('ls', [tcommons.data_folder.str()], ['ls']),
            (['ls'], [tcommons.data_folder.str()], ['ls']),
            (['ls -ltr'], [tcommons.data_folder.str()], ['ls', '-ltr']),
            ('ls -ltr', [tcommons.data_folder.str()], ['ls', '-ltr']),
            (['ls', '-ltr'], [tcommons.data_folder.str()], ['ls', '-ltr']),
            ],
        )
    def test_SubprocessTask_1(self, in1, in2, expected):
        """Test cmd_exec from string."""
        sub = LM.SubprocessTask(in1,  in2)
        assert sub.cmd_exec == expected

    def test_SubprocessTask_2(self):
        """Test prepare_cmd()."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        sub.prepare_cmd()
        assert sub.cmd == ['ls', tcommons.data_folder.str()]


    def test_SubprocessTask_3(self):
        """Test execute."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        sub.prepare_cmd()
        sub.execute()


    def test_SubprocessTask_4(self):
        """Test result is CompletedProcess."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        sub()
        assert isinstance(sub.result, subprocess.CompletedProcess)


    def test_SubprocessTask_5(self):
        """Test string before prepare_cmd()."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        assert str(sub) == "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
            tcommons.data_folder.str())


    def test_SubprocessTask_6(self):
        """Test repr()."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        assert repr(sub) == \
            "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
                tcommons.data_folder.str())


    def test_SubprocessTask_7(self):
        """Test str() after prepare_cmd()."""
        sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
        sub.prepare_cmd()
        assert str(sub) == 'ls {}'.format(tcommons.data_folder.str())


    def test_SubprocessTask_8(self):
        """Test input None."""
        sub = LM.SubprocessTask('ls')
        sub()


class TestDSSPTask:

    @pytest.mark.parametrize(
        'in1,in2',
        [
            ('dssp', '1XXX.pdb'),
            (['dssp'], '1XXX.pdb'),
            ],
        )
    def test_DSSPTask_1(self, in1, in2):
        """Test DSSPTask init."""
        LM.DSSPTask(in1, in2)

    def test_DSSPTask_2(self):
        """Test raises ValueError when input is missing."""
        with pytest.raises(TypeError):
            LM.DSSPTask('dssp')
  
    @pytest.mark.parametrize(
        'in1,in2,expected',
        [
            ('dssp', '1XXX.pdb', ['dssp', '-i', '1XXX.pdb']),
            (['dssp'], '1XXX.pdb', ['dssp', '-i', '1XXX.pdb']),
            ]
        )
    def test_DSSPTask_3(self, in1, in2, expected):
        """Test DSSPTask prepare_cmd."""
        dssptask = LM.DSSPTask(in1, in2)
        dssptask.prepare_cmd()
        assert dssptask.cmd == expected

    def test_DSSPTask_pdb_path(self):
        dssptask = LM.DSSPTask('dssp', '1XXX.pdb')
        assert dssptask.pdb_path == '1XXX.pdb'

    def test_DSSPTask_call(self):
        dssptask = LM.DSSPTask('ls', '-ltr')
        dssptask()


class TestWorker:

    def test_worker_1(self):
        LM.Worker(None, None)

    def test_worker_2(self):
        with pytest.raises(TypeError):
            LM.Worker()

    def test_worker_3(self):
        """Test Worker run."""
        tasks = multiprocessing.JoinableQueue()
        queue = multiprocessing.Queue()
        workers = [LM.Worker(tasks, queue, timeout=1)]
        for w in workers:
            w.start()
        tasks.put(LM.SubprocessTask('ls -ltr'))
        tasks.put(None)
        tasks.join()
        numjobs = 1
        while numjobs:
            queue.get()
            numjobs -= 1


class TestJoinedResults:
    
    def test_init(self):
        jr = LM.JoinedResults(
            '-ltr',
            'ls',
            ncores=7,
            TaskMethod=LM.Task,
            results_parser=int)
        assert jr.input_data == '-ltr'
        assert jr.cmd == 'ls'
        assert jr.ncores == 7
        assert jr.TaskMethod == LM.Task
        assert jr.rp == int

    def test_build(self):
        jr = LM.JoinedResults('-ltr', 'ls', ncores=7)
        jr.build()
        assert isinstance(jr.tasks, multiprocessing.queues.JoinableQueue)
        assert isinstance(jr.queue_results, multiprocessing.queues.Queue)
        assert len(jr.workers) == 7
        assert all(isinstance(w, LM.Worker) for w in jr.workers)

    def test_run(self):
        jr = LM.JoinedResults('-ltr', 'ls', ncores=7)
        jr.run()
