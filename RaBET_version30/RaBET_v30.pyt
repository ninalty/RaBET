import os
import sys

sys.dont_write_bytecode = True
scripts_dir = os.path.join(os.path.dirname(__file__), 'scripts')
sys.path.append(scripts_dir)

from WCMapTool import WCMapTool

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "RaBET"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [WCMapTool]
