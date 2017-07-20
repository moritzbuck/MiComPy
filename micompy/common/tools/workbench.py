from micompy.common.tools.bbmap import BBmap
from micompy.common.tools.checkm import CheckM
from micompy.common.tools.mash import MASH
from micompy.common.tools.hmmer import HMMer
from micompy.common.tools.tool import Tool

class WorkBench(object):
    def __getitem__(self, key):
        return self.tools.get(key)

    def __init__(self):
        self.tools = {}

    def add_tool(self, tool, name = None):
        assert isinstance(tool, Tool), tool + " is not a tool"
        self.tools[name if name else tool.name] = tool

    def default_bench(self):
        self.add_tool(BBmap())
        self.add_tool(CheckM())
        self.add_tool(MASH())
        self.add_tool(HMMer())
