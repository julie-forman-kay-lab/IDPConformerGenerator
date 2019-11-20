"""Multi-core related objects."""
import multiprocessing
import subprocess

from idpconfgen import Path, log
from idpconfgen.logger import S, T
from idpconfgen.libs import libcheck


class Worker(multiprocessing.Process):
    
    def __init__(self, task_queue, result_queue, timeout=10):
        
        # super().__init__(self) ?
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.timeout = timeout

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                log.info(S('exiting: {}', proc_name))
                self.task_queue.task_done()
                break
            result = next_task()
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class Task:
    pass


class SubprocessTask(Task):
    
    @libcheck.argstype(Task, (list, str))
    @libcheck.kwargstype((type(None), list, tuple))
    def __init__(self, cmd, input_=None):
        """
        General task operation.

        Parameters
        ----------
        cmd : str or list
            The command to execute. Options to the command
            should be given here, and two formats are possible,
            as string or as list. If string, string is split into
            its components, for example:
                
                'ls -ltr' will result in ['ls', '-ltr']
            
            A preprepared list can be given instead.
        
        input : list
            A list containing whatever input the command must receive.
            Defaults to None, no input is used.
        """
        try:
            self.cmd_exec = cmd.split()
        except AttributeError:  # in case cmd is a list
            self.cmd_exec = cmd
        
        self.input = input_


    def __str__(self):
        try:
            return ' '.join(self.cmd)
        except AttributeError:
            return repr(self)

    def __repr__(self):
        return '{}({})'.format(
            __class__.__name__,
            ','.join('{}={}'.format(k, v) for k, v in self.__dict__.items()),
            )
    
    def __call__(self):
        self.prepare_cmd()
        self.execute()
        return self.result.stdout.decode('utf8')

    def prepare_cmd(self):
        try:
            self.cmd = self.cmd_exec + self.input
        except TypeError:  # in case input_ is None
            self.cmd = self.cmd_exec

    def execute(self):
        log.info(S('running {}', self.cmd))
        self.result = subprocess.run(
            self.cmd,
            capture_output=True,
            )


class DSSPTask(SubprocessTask):
    """Subprocess Task for DSSP third party executable."""
    
    @libcheck.argstype(Task, (list, str), str)
    def __init__(self, cmd, input_):
        # forces input_ to be positional parameter in subclass DSSPTask
        
        try:
            cmd.append('-i')
        except AttributeError:
            cmd = [cmd, '-i']

        self.pdb_path = input_
        input_ = [input_]

        super().__init__(cmd, input_=input_)
    
    def __call__(self):
        self.prepare_cmd()
        self.execute()
        return (self.pdb_path, self.result.stdout.decode('utf8'))


class JoinedResults:
    
    def __init__(
            self,
            input_data,
            command,
            ncores=1,
            TaskMethod=SubprocessTask,
            results_parser=str,
            ):
        self.input_data = input_data
        self.cmd = command
        assert isinstance(self.cmd, str), '{}'.format(type(self.cmd))
        self.ncores = 1
        self.TaskMethod = TaskMethod
        self.rp = results_parser
        
    def build(self):
        """Build tasks, results and workers."""
        self.tasks = multiprocessing.JoinableQueue()
        self.queue_results = multiprocessing.Queue()
        self.workers = \
            [Worker(self.tasks, self.queue_results)
                for i in range(self.ncores)]
    
    def run(self):
        """Run."""
        self.build()
        
        for w in self.workers:
            w.start()

        for input_datum in self.input_data:
            self.tasks.put(
                self.TaskMethod(self.cmd, input_datum))

        # Adds a poison pill
        for w in self.workers:
            self.tasks.put(None)

        self.tasks.join()

        numjobs = len(self.input_data)
        self.results = []
        while numjobs:
            self.results.append(self.rp(self.queue_results.get()))
            numjobs -= 1
        print(self.results)
