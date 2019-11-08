"""Test libs for multicore operations."""
import pytest
import subprocess

from idpconfgen.libs import libmulticore as LM

from . import tcommons

def test_Task_1():
    LM.Task()


def test_SubprocessTask_1():
    """Test cmd_exec from string."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert sub.cmd_exec == ['ls']


def test_SubprocessTask_1_1():
    """Test cmd_exec from list."""
    sub = LM.SubprocessTask(['ls'],  [tcommons.data_folder.str()])
    assert sub.cmd_exec == ['ls']


def test_SubprocessTask_1_2():
    """Test cmd_exec from list with parameters."""
    sub = LM.SubprocessTask(['ls', '-ltr'],  [tcommons.data_folder.str()])
    assert sub.cmd_exec == ['ls', '-ltr']


def test_SubprocessTask_1_3():
    """Test cmd_exec from string with parameters."""
    sub = LM.SubprocessTask('ls -ltr',  [tcommons.data_folder.str()])
    assert sub.cmd_exec == ['ls', '-ltr']


def test_SubprocessTask_2():
    """Test prepare_cmd()."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    assert sub.cmd == ['ls', tcommons.data_folder.str()]


def test_SubprocessTask_3():
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    sub.execute()


def test_SubprocessTask_4():
    """Test result is CompletedProcess."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub()
    assert isinstance(sub.result, subprocess.CompletedProcess) 


def test_SubprocessTask_5():
    """Test string before prepare_cmd()."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert str(sub) == "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
        tcommons.data_folder.str())


def test_SubprocessTask_6():
    """Test repr()."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    assert repr(sub) == "SubprocessTask(cmd_exec=['ls'],input=['{}'])".format(
        tcommons.data_folder.str())


def test_SubprocessTask_7():
    """Test str() after prepare_cmd()."""
    sub = LM.SubprocessTask('ls',  [tcommons.data_folder.str()])
    sub.prepare_cmd()
    assert str(sub) == 'ls {}'.format(tcommons.data_folder.str())


def test_SubprocessTask_8():
    """Test input None."""
    sub = LM.SubprocessTask('ls')
    sub()


def test_DSSPTask_1():
    """Test DSSPTask init."""
    LM.DSSPTask('dssp', '1XXX.pbd')


def test_DSSPTask_1_1():
    """Test cmd list."""
    LM.DSSPTask(['dssp'], '1XXX.pbd')


def test_DSSPTask_2():
    """Test raises ValueError when input is missing."""
    with pytest.raises(TypeError):
        LM.DSSPTask('dssp')
        

def test_DSSPTask_3():
    """Test DSSPTask prepare_cmd."""
    dssptask = LM.DSSPTask('dssp', '1XXX.pdb')
    dssptask.prepare_cmd()
    assert dssptask.cmd == ['dssp', '-i', '1XXX.pdb']


def test_DSSPTask_4():
    """Test DSSPTask prepare_cmd from list."""
    dssptask = LM.DSSPTask(['dssp'], '1XXX.pdb')
    dssptask.prepare_cmd()
    assert dssptask.cmd == ['dssp', '-i', '1XXX.pdb']
