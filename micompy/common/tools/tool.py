import subprocess
import os

class Tool(object):

    def __init__(self, name, executable):
        self.name = name
        self.executable = executable

        try:
            with open(os.devnull, 'w') as handle:
                subprocess.call([self.executable, "-h"], stdout=handle)
        except OSError as e:
            print("%s's executable %s not found" % (self.name, self.executable))
