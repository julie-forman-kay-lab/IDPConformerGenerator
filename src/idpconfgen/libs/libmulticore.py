"""Multi-core related objects."""
import subprocess

from idpconfgen.logger import S, T
from idpconfgen.libs import libcheck


class Worker:
    
    def __init___(self, task_queue, result_queue, timeout=10):
        
        # super().__init__(self) ?
        multiprocessing.Process.__init__(self)
        self.task_queque = task_queue
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
            result = next_tasks()
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class Task:
    pass


class SubprocessTask(Task):
    
    @libcheck.argstype(Task, (list, str))
    @libcheck.kwargstype((type(None), list, tuple),)
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
            cmd_exec = cmd.split()
        except AttributeError:  # in case cmd is a list
            cmd_exec = cmd

        try:
            self.cmd = cmd_exec.extend(input_)
        except TypeError:  # in case input_ is None
            pass

        def __str__(self):
            return ' '.join(self.cmd)

        def __call__(self):
            log.info(S('running', self.cmd))
            result = subprocess.run(
                self.cmd,
                capture_output=True,
                )
            return result.stdout  #.decode('utf8')


class JoinedResults:
    
    def __init__(
            self,
            input_data,
            command,
            ncores=1,
            task_method=SubprocessTask,
            results_parser=str,
            ):
        self.input_data = input_data
        self.cmd = command
        self.ncores = 1
        self.task_method = taskmethod
        self.rp = results_parser
        
    def build(self):
        """Build tasks, results and workers."""
        self.tasks = multiprocessing.JoinableQueue()
        self.results = multiprocessing.Queue()
        self.workers = [Worker(tasks, results) for i in range(self.ncores)]
    
    def run(self):
        """Run."""
        self.build()
        
        for w in self.workers:
            w.start()

        for input_datum in self.input_data:
            self.tasks.put(CmdTask(self.cmd, input_datum))

        # Adds a poison pill
        for w in self.workers:
            tasks.put(None)

        tasks.join()

        numjobs = len(self.input_data)
        self.results = []
        while numjobs:
            self.results.append(self.rp(self.results.get()))
            numjobs -= 1
