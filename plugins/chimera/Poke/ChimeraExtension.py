import chimera.extension

class Poke(chimera.extension.EMO):
	def name(self):
		return 'Poke scouter'
	def description(self):
		return 'Find pocket on a protein'
	def categories(self):
		return ['Structure Analysis']
	def activate(self):
		from chimera.dialogs import display
		display(self.module('gui').PokeDialog.name)
		return None
	def modelPanelCB(self, molecules):
		from chimera.dialogs import display
		d = display(self.module('gui').PokeDialog.name)
		d.molList.setvalue(molecules)

emo = Poke(__file__)
chimera.extension.manager.registerExtension(emo)
import ModelPanel
ModelPanel.addButton("find pockets...", emo.modelPanelCB,
						defaultFrequent=False)

def cmdPoke(cmdName, args):
    #     from Midas.midas_text import doExtensionFunc
	# import AddH
    #     doExtensionFunc(AddH.cmdAddH, args,
    #                             specInfo=[("spec", "molecules", "molecules")])
	pass

from Midas.midas_text import addCommand
addCommand("poke", cmdPoke, help=True)

