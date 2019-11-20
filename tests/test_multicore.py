"""Test libs for multicore operations."""
import pytest
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
