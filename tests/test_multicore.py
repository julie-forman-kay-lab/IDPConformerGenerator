"""Test libs for multicore operations."""
import subprocess

from idpconfgen.libs import libmulticore as LM

from . import tcommons

def test_Task_1():
    LM.Task()


def test_SubprocessTask_1():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert sub.cmd_exec == ['ls']

def test_SubprocessTask_2():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    assert sub.cmd == ['ls', tcommons.data_folder.str()]

def test_SubprocessTask_3():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    sub.execute()

def test_SubprocessTask_4():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub()
    assert isinstance(sub.result, subprocess.CompletedProcess) 

def test_SubprocessTask_5():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert str(sub) == "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
        tcommons.data_folder.str())

def test_SubprocessTask_6():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert repr(sub) == "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
        tcommons.data_folder.str())

def test_SubprocessTask_7():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    assert str(sub) == 'ls {}'.format(tcommons.data_folder.str())
