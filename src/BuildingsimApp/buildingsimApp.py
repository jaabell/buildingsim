#!/usr/bin/python

from Tkinter import *
import tkSimpleDialog
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfilename
#import elementtree.ElementTree as ET
import sys
import os
python_version = sys.version_info[:2]
if python_version >= (2, 5):
	import xml.etree.cElementTree as ET
else:
	try:
		import cElementTree as ET
	except ImportError:
		try:
			import ElementTree as ET
		except ImportError:
			msg = "\nYou need the [c]ElementTree package\n" \
			"from http://effbot.org/downloads/\n\n";
			sys.stderr.write(msg)
			raise
			print >> sys.stderr, "ET imported from", ET.__file_
from buildingtool import *

class App(Frame):
    def run_action(self):
        building = {}
        building['dx'] = float(self.dxEntry.get())     	# [m]
        building['dy'] =  float(self.dyEntry.get())     # [m]
        building['g'] = 9.806				# Gravitational Constant [m/s^2]
        building['xsi'] = 5. 				# Modal damping  [%]

        building['h'] = zeros((len(self.buildings,)))       # Altura entrepiso [m]
        building['E'] = zeros((len(self.buildings),))       # Modulo elasticidad piso [tonf/cm**2]
        building['I'] = zeros((len(self.buildings),))       # Momento inercia columnas [m**4]
        building['gamma'] = zeros((len(self.buildings),))   # Peso unitario losas [tonf/m**3]
        building['esp'] = zeros((len(self.buildings),))     # espesor losas [m]

        for i, item in enumerate(self.buildings):
            for key,value in item.iteritems():
                building[key][i] = float(value)
        
        building['name'] = self.buildingName

        building_sim(building)
        #tkMessageBox.showinfo("Procesando", "Estamos procesando por usted...")

    def add_floor(self):
        floor = {"h":self.hEntry.get(), "E":self.eEntry.get(), "I":self.iEntry.get(), "gamma":self.gammaEntry.get(), "esp":self.espEntry.get()}
        self.buildings.append(floor)
        self.update_floor_listbox()

    def edit_floor(self):
        items = map(int, self.floorListBox.curselection())
        if len(items) > 0:
            b = self.buildings[items[0]]
            self.hEntry.delete(0, END)
            self.hEntry.insert(0, b["h"])
            self.eEntry.delete(0, END)
            self.eEntry.insert(0, b["E"])
            self.iEntry.delete(0, END)
            self.iEntry.insert(0, b["I"])
            self.gammaEntry.delete(0, END)
            self.gammaEntry.insert(0, b["gamma"])
            self.espEntry.delete(0, END)
            self.espEntry.insert(0, b["esp"])

    def change_floor(self):
        items = map(int, self.floorListBox.curselection())
        if len(items) > 0:
            self.buildings[items[0]] = {"h":self.hEntry.get(), "E":self.eEntry.get(), "I":self.iEntry.get(), "gamma":self.gammaEntry.get(), "esp":self.espEntry.get()}
        self.update_floor_listbox()

    def delete_floor(self):
        items = map(int, self.floorListBox.curselection())
        for floor in items:
            self.buildings.remove(self.buildings[floor])
        self.update_floor_listbox()

    def update_floor_listbox(self):
        self.floorListBox.delete(0, END)
        for item in self.buildings:
            s = ""
            for key,value in item.iteritems():
                s = s + "[" + key + "=" + value + "]"
            self.floorListBox.insert(END, s)

    def prepare_to_building(self):
        if self.buildingName:
            self.buildings = []
            Label(self, text="Edificio ", font=("Verdana", 14), bg=self.mainBg).grid(row=0, column=0, padx=5, pady=5, sticky=NE)
            self.nameLabel = Label(self, text=self.buildingName, bg=self.mainBg)
            self.nameLabel["fg"] = "red"
            self.nameLabel["font"] = ("Verdana", 14)
            self.nameLabel.grid(row=0, column=1, padx=5, pady=5, sticky=NW)

            # Accion
            self.runButton = Button(self, text="Calcular", command=self.run_action, font=("Verdana", 14), width=16, height=2)
            self.runButton.grid(row=0, column=3, columnspan=2, padx=10, pady=5)

            # Base
            Label(self, text="Base", font=("Verdana", 12), bg=self.mainBg).grid(row=1, column=0, columnspan=3, padx=5, pady=10)
            Label(self, text="Largo [m]", font=("Verdana", 10), bg=self.mainBg).grid(row=2, column=0, padx=5, sticky=E)
            Label(self, text="Ancho [m]", font=("Verdana", 10), bg=self.mainBg).grid(row=3, column=0, padx=5, sticky=E)

            self.dxEntry = Entry(self)
            self.dxEntry.grid(row=2, column=1, padx=3)

            self.dyEntry = Entry(self)
            self.dyEntry.grid(row=3, column=1, padx=3)

            # Pisos
            Label(self, text="Pisos", font=("Verdana", 12), bg=self.mainBg).grid(row=4, column=0, columnspan=3, padx=5, pady=5)
            Label(self, text="Altura [m]", font=("Verdana", 10), bg=self.mainBg).grid(row=5, column=0, padx=5, sticky=E)
            self.hEntry = Entry(self)
            self.hEntry.grid(row=5, column=1, padx=5)

            Label(self, text="Modulo Elasticidad [tonf/m^2]", font=("Verdana", 10), bg=self.mainBg).grid(row=6, column=0, padx=5, sticky=E)
            self.eEntry = Entry(self)
            self.eEntry.grid(row=6, column=1, padx=5)

            Label(self, text="Momento Inercia [m^4]", font=("Verdana", 10), bg=self.mainBg).grid(row=7, column=0, padx=5, sticky=E)
            self.iEntry = Entry(self)
            self.iEntry.grid(row=7, column=1, padx=5)

            Label(self, text="Peso Unitario Losa [tonf/m^3]", font=("Verdana", 10), bg=self.mainBg).grid(row=8, column=0, padx=5, sticky=E)
            self.gammaEntry = Entry(self)
            self.gammaEntry.grid(row=8, column=1, padx=5)

            Label(self, text="Espesor [m]", font=("Verdana", 10), bg=self.mainBg).grid(row=9, column=0, padx=5, sticky=E)
            self.espEntry = Entry(self)
            self.espEntry.grid(row=9, column=1, padx=5)

            self.addFloorButton = Button(self, text="Agregar", command=self.add_floor, width=10)
            self.addFloorButton.grid(row=5, column=2, rowspan=2, padx=5)

            self.editFloorButton = Button(self, text="Cambiar", command=self.change_floor, width=10)
            self.editFloorButton.grid(row=6, column=2, rowspan=2, padx=5)

            self.floorListBox = Listbox(self, width=40)
            self.floorListBox.grid(row=5, column=3, rowspan=4, columnspan=2, padx=10)

            self.deleteFloorButton = Button(self, text="Editar", command=self.edit_floor, width=10)
            self.deleteFloorButton.grid(row=9, column=3, padx=5)

            self.deleteFloorButton = Button(self, text="Eliminar", command=self.delete_floor, width=10)
            self.deleteFloorButton.grid(row=9, column=4, padx=5)

    def prepare_new_building(self):
        self.buildingName = tkSimpleDialog.askstring("Ingreso", "Ingresa el nombre del Edificio:", parent=self)
        self.prepare_to_building()

    def load_building(self):
        filename = askopenfilename(title="Elige un archivo para abrir", initialdir=".", filetypes=[("Buildingsim File","*.bsim")])
        tree = ET.parse(filename)
        root = tree.getroot()
        self.buildingName = root.get("name")        
        self.prepare_to_building()
        
        self.dxEntry.insert(0, root.get("dx"))
        self.dyEntry.insert(0, root.get("dy"))
        
        for floorNode in tree.getiterator("floor"):
            floor = {"h":floorNode.get("h"), "E":floorNode.get("E"), "I":floorNode.get("I"), "gamma":floorNode.get("gamma"), "esp":floorNode.get("esp")}
            self.buildings.append(floor)
            self.update_floor_listbox()
            
    def save_building(self):
        filename = asksaveasfilename(parent=self.master, title="Ingresa o selecciona el archivo donde guardar", filetypes=[("Buildingsim File","*.bsim")])

        buildingNode = ET.Element("building")
        buildingNode.set("name", self.buildingName)
        buildingNode.set("dx", self.dxEntry.get())
        buildingNode.set("dy", self.dyEntry.get())

        floorsListNode = ET.SubElement(buildingNode, "floors")

        #Save floors
        for i, item in enumerate(self.buildings):
            floorNode = ET.SubElement(floorsListNode, "floor")
            for key,value in item.iteritems():
                floorNode.set(key, value)

        tree = ET.ElementTree(buildingNode)
        tree.write(filename)

    def create_menu(self):
        self.menu = Menu(self)
        self.master.config(menu=self.menu)

        self.buildingMenu = Menu(self.menu)
        self.buildingMenu.add_command(label="Nuevo", command=self.prepare_new_building)
        self.buildingMenu.add_command(label="Abrir", command=self.load_building)
        self.buildingMenu.add_separator()
        self.buildingMenu.add_command(label="Guardar", command=self.save_building)
        self.buildingMenu.add_separator()
        self.buildingMenu.add_command(label="Salir", command=self.quit)
        self.menu.add_cascade(label="Edificio", menu=self.buildingMenu)

    def __init__(self, master):
        self.mainBg = "white"
        Frame.__init__(self, master, width=master.winfo_screenwidth(), height=master.winfo_screenheight(), bg=self.mainBg)
        master.config(bg="#003C6E")
        self.create_menu()
        
        self.grid(padx=20, pady=20,sticky=N+S+E+W)
        self.grid_propagate(0)
        
        self.master.title("Calculadora de Comportamiento de Edificios")        
        
        self.canvas = Canvas(width = 77, height = 98)        
        self.gifLogo = PhotoImage(file = 'imgs/ing.gif')
        self.canvas.create_image(0, 0, image = self.gifLogo, anchor = NW)
        self.canvas.grid(row=0, column=0, padx=0, pady=0, sticky=NW)

#main
root = Tk()
toplevel = root.winfo_toplevel()
# Windows
try:    
    toplevel.wm_state('zoomed')
#Others
except:
    w, h = root.winfo_screenwidth(), root.winfo_screenheight() - 60
    geom_string = "%dx%d+0+0" % (w,h)
    toplevel.wm_geometry(geom_string)
exeFileName = sys.argv[0]
exeDirectory = os.path.dirname(exeFileName)
os.chdir(exeDirectory)
w = App(root)
root.mainloop()
