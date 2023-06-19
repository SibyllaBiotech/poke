import chimera
from chimera.baseDialog import BaseDialog
from chimera.widgets import MoleculeScrolledListBox
import Tkinter
import Tix

class PokeDialog(BaseDialog):
	name = "Pocket scouter"
	buttons = ('Close',)
	default = 'Close'

	def __init__(self, **kw):
		self.parent = None

		self.scrolledHListPockets = None
		self.hlistPockets = None

		self.pockets = []
		self.descriptors_names = []

		chimera.triggers.addHandler('Atom', 
					 self.checkExistingPockets, None)

		BaseDialog.__init__(self, **kw)

	def fillInUI(self, parent):
		self.parent = parent
		parent.grid_columnconfigure(0, minsize=1000, weight=1000)

		row = 0		

		l = Tkinter.Label(parent, text=
		    'Select the protein to search for pockets on.\n'
		    'Be sure that the protein is protonated before continuing.')
		l.grid(row=row, column=0, sticky="news", pady=(30, 0), padx=(10, 10))
		row += 1

		self.runButton = Tkinter.Button(parent,
					text='Find pockets on selection' ,
					command=self.run_poke)
		self.runButton.grid(row=row, column=0, sticky="news", padx=(10, 10), pady=(20, 20))
		row += 1

		l = Tkinter.Label(parent, text="Pockets")
		l.configure(font=('Helvatical bold',12))
		l.grid(row=row, column=0, sticky="news", pady=(50, 10), padx=(10, 10))
		row += 1

		self.updateDisplayedPockets()

	def run_poke(self):
		import os
		import subprocess
		import shutil
		import csv

		try:
			with open('chimera_poke.pdb', 'w') as file:
				serial = 1
				for chains in chimera.selection.currentChains():
					for residue in chains.residues:
						for atom in residue.atoms:
							file.write(("ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
										"{chainID:1s}{resSeq:4d}{iCode:1s}"
										"   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}"
									#	"{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"
										"\n").format(
								serial=serial, name=atom.name, altLoc=atom.altLoc, 
								resName=atom.residue.type, chainID='A', resSeq=atom.residue.id.position,
								iCode=' ',
								x=atom.coord().data()[0],
								y=atom.coord().data()[1],
								z=atom.coord().data()[2],
								occupancy=0))
							serial += 1

				if serial == 1:
					raise RuntimeError('The selection is empty.')


			try:
				subprocess.check_output(['poke', 'chimera_poke.pdb', '-dv'], stderr=subprocess.STDOUT)
			except subprocess.CalledProcessError as error:
				raise RuntimeError(error.output)
			
			mList = chimera.openModels.open('chimera_poke/volumes.pdb',
						identifyAs='Pockets')

			chimera.runCommand('rangecolor bfactor -1 blue 0 white +1 red #'+ str(mList[0].id))

			with open('chimera_poke/descriptors.csv', 'rb') as csvfile:
				reader = csv.reader(csvfile)
				for row, values in enumerate(reader):
					if row == 0:
						self.descriptors_names = list(values)[2:-2]
					else:
						self.pockets.append(Pocket(
							model_id=mList[0].id,
							resid=800 + row - 1,
							descriptors=values[2:-2]))
						
			if row == 0:
				chimera.openModels.close(mList)

			self.updateDisplayedPockets()
											
		finally:
			try:
				os.unlink('chimera_poke.pdb')
			except:
				pass
			try:
				shutil.rmtree('chimera_poke', ignore_errors=True)
			except:
				pass

	def updateDisplayedPockets(self):
		self.scrolledHListPockets = Tix.ScrolledHList(self.parent,
				options = """hlist.columns %d
				hlist.header 1
				hlist.selectMode extended
				hlist.indicator 1""" % len(self.descriptors_names))
		self.hlistPockets = self.scrolledHListPockets.hlist
		self.hlistPockets.config(browsecmd=self.updateSelectedPockets)

		style = Tix.DisplayStyle('text', refwindow=self.scrolledHListPockets, anchor='c')

		for cix, cname in enumerate(self.descriptors_names):
			self.hlistPockets.header_create(cix, itemtype='window',
								window=ColumnHeader(self.hlistPockets, ' ' + cname + ' ', command=lambda e: 0))

		for rix, pock in enumerate(self.pockets):
			self.hlistPockets.add(rix)
			for cix, cv in enumerate(pock.descriptors):
				self.hlistPockets.item_create(rix, cix, text=cv, style=style)

		self.scrolledHListPockets.grid(row=4, column=0, sticky="news", padx=(10, 10))

	def updateSelectedPockets(self, event=None):
		selection = ','.join([self.pockets[int(s)].selector() for s in self.hlistPockets.info_selection()])
		chimera.runCommand('select '+selection)

	def checkExistingPockets(self, trigName, myData, atoms):
		if len(atoms.deleted):
			existing_pockets = []
			for p in self.pockets:
				sel = chimera.selection.OSLSelection(p.selector()).atoms()	
				if len(sel) > 0:
					existing_pockets.append(p)

			if len(self.pockets) != len(existing_pockets):
				self.pockets = existing_pockets
				self.updateDisplayedPockets()



class Pocket:
	def __init__(self, model_id, resid, descriptors):
		self.model_id = model_id
		self.resid = resid
		self.descriptors = descriptors

	def selector(self):
		return '#' + str(self.model_id) + ':' + str(self.resid)

class ColumnHeader(Tkinter.Frame):

	def __init__(self, master, title, **kw):
		try:
			self.colCmd = kw["command"]
		except KeyError:
			self.colCmd = None
		else:
			del kw["command"]
		kw["takefocus"] = 1
		Tkinter.Frame.__init__(self, master)
		self.colTitle = Tkinter.Label(self, text=title, **kw)
		self.colTitle.pack(side=Tkinter.LEFT)
		self.colImage = None
		if self.colCmd:
			self.colTitle.bind("<ButtonRelease>", self.colCmd)

	def setTitle(self, title):
		self.colTitle.config(text=title)

	def setImage(self, img):
		if img is None:
			if self.colImage:
				self.colImage.destroy()
				self.colImage = None
		else:
			if self.colImage:
				self.colImage.config(image=img)
			else:
				self.colImage = Tkinter.Label(self, image=img)
				self.colImage.pack(side=Tkinter.RIGHT)
				if self.colCmd:
					self.colImage.bind("<ButtonRelease>",
								self.colCmd)

from chimera import dialogs
dialogs.register(PokeDialog.name, PokeDialog)
