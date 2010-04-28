#!/usr/bin/python

from Tkinter import *
import tkSimpleDialog
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfilename
from ttk import *
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
from definitions import *

class App(PanedWindow):
    def prepare_action(self):
        self.building = {}
        self.building['dx'] = float(self.dxEntry.get())     	# [m]
        self.building['dy'] =  float(self.dyEntry.get())     # [m]
        self.building['g'] = 9.806				# Gravitational Constant [m/s^2]
        self.building['xsi'] = 5. 				# Modal damping  [%]

        self.building['h'] = zeros((len(self.buildings,)))       # Altura entrepiso [m]
        self.building['E'] = zeros((len(self.buildings),))       # Modulo elasticidad piso [tonf/cm**2]
        self.building['I'] = zeros((len(self.buildings),))       # Momento inercia columnas [m**4]
        self.building['gamma'] = zeros((len(self.buildings),))   # Peso unitario losas [tonf/m**3]
        self.building['esp'] = zeros((len(self.buildings),))     # espesor losas [m]

        for i, item in enumerate(self.buildings):
            for key,value in item.iteritems():
                self.building[key][i] = float(value)
        
        self.building['name'] = self.buildingName

        self.building = form(self.building)
        
        #Enable sim and plot tabs
        self.actionsTabs.tab(1, state="normal")
        self.actionsTabs.tab(2, state="normal")        
        self.actionsTabs.tab(3, state="normal")
    
    def disabled_tabs(self):
        self.actionsTabs.tab(1, state="disabled")
        self.actionsTabs.tab(2, state="disabled")        
        self.actionsTabs.tab(3, state="disabled")
        self.actionsTabs.tab(4, state="disabled")
        self.actionsTabs.tab(5, state="disabled")
        self.actionsTabs.tab(6, state="disabled")
    
    def add_floor(self):
        floor = {"h":self.hEntry.get(), "E":self.eEntry.get(), "I":self.iEntry.get(), "gamma":self.gammaEntry.get(), "esp":self.espEntry.get()}
        self.buildings.append(floor)
        self.update_floor_listbox()
        #Disable sim and plot tabs
        self.disabled_tabs()

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
            #Disable sim and plot tabs
            self.disabled_tabs()
        self.update_floor_listbox()

    def delete_floor(self):
        items = map(int, self.floorListBox.curselection())
        for floor in items:
            self.buildings.remove(self.buildings[floor])
        self.update_floor_listbox()
        #Disable sim and plot tabs
        self.disabled_tabs()

    def update_floor_listbox(self):
        self.floorListBox.delete(0, END)
        for item in self.buildings:
            s = ""
            for key,value in item.iteritems():
                s = s + "[" + key + "=" + value + "]"
            self.floorListBox.insert(END, s)

    def prepare_to_building(self):
        self.buildings = []
        if len(self.buildingName) != 0 and self.isCreated == False:
            self.isCreated = True            
            self.create_building_frame()
            self.create_actions_frame()
    
    def create_building_frame(self):
        currentRow = 0
        #Create Widgets            
        Label(self.frameLeft, text="Edificio ", style="Title.TLabel").grid(row=currentRow, column=0, padx=5, pady=5, sticky=NE)
        self.nameLabel = Label(self.frameLeft, text=self.buildingName, style="Red.TLabel")
        self.nameLabel.grid(row=currentRow, column=1, padx=5, pady=5, sticky=NW)

        # Base
        currentRow += 1
        Label(self.frameLeft, text="Base", style="Title.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=10)
        currentRow += 1
        Label(self.frameLeft, text="Largo [m]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.dxEntry = Entry(self.frameLeft)
        self.dxEntry.grid(row=currentRow, column=1, padx=3)
        currentRow += 1
        Label(self.frameLeft, text="Ancho [m]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.dyEntry = Entry(self.frameLeft)
        self.dyEntry.grid(row=3, column=1, padx=3)
        
        currentRow += 1

        # Pisos
        Label(self.frameLeft, text="Pisos", style="Title.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=5)
        currentRow += 1
        Label(self.frameLeft, text="Altura [m]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.hEntry = Entry(self.frameLeft)
        self.hEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        Label(self.frameLeft, text="Modulo Elasticidad [tonf/m^2]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.eEntry = Entry(self.frameLeft)
        self.eEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        Label(self.frameLeft, text="Momento Inercia [m^4]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.iEntry = Entry(self.frameLeft)
        self.iEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        Label(self.frameLeft, text="Peso Unitario Losa [tonf/m^3]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.gammaEntry = Entry(self.frameLeft)
        self.gammaEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        Label(self.frameLeft, text="Espesor [m]").grid(row=currentRow, column=0, padx=5, sticky=E)
        self.espEntry = Entry(self.frameLeft)
        self.espEntry.grid(row=currentRow, column=1, padx=5)
                    
        self.addFloorButton = Button(self.frameLeft, text="Agregar >>", command=self.add_floor, width=10)
        self.addFloorButton.grid(row=6, column=2, rowspan=2, padx=5)

        self.editFloorButton = Button(self.frameLeft, text="Cambiar", command=self.change_floor, width=10)
        self.editFloorButton.grid(row=7, column=2, rowspan=2, padx=5)

        self.floorListBox = Listbox(self.frameLeft, width=40)
        self.floorListBox.grid(row=6, column=3, rowspan=4, columnspan=2, padx=10)

        self.deleteFloorButton = Button(self.frameLeft, text="Editar", command=self.edit_floor, width=10)
        self.deleteFloorButton.grid(row=10, column=3, padx=5)

        self.deleteFloorButton = Button(self.frameLeft, text="Eliminar", command=self.delete_floor, width=10)
        self.deleteFloorButton.grid(row=10, column=4, padx=5)

    def create_actions_frame(self):
        self.actionsTabs = Notebook(self.frameRight)
        self.actionsTabs.pack(fill=BOTH, expand=1)
        self.create_prepare_action_frame()
        self.create_sim_action_frame()
        self.create_freqresp_action_frame()
        self.create_compare_action_frame()
        self.create_envresp_action_frame()
        self.create_floorresp_action_frame()
        self.create_animation_action_frame()
    
    def create_prepare_action_frame(self):
        self.prepareFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.prepareFrame, text="Preparar")
        
        # Prepare Button
        self.prepareButton = Button(self.prepareFrame, text="Preparar", command=self.prepare_action)
        self.prepareButton.grid(row=0, column=0, padx=10, pady=10)
        
    def create_sim_action_frame(self):
        self.simFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.simFrame, text="Simular", state="disabled")
        
    def create_freqresp_action_frame(self):
        self.freqrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.freqrespFrame, text="Resp. a Frecuencia", state="disabled")
        
    def create_compare_action_frame(self):
        self.compareFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.compareFrame, text="Comparar", state="disabled")
    
    def create_envresp_action_frame(self):
        self.envrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.envrespFrame, text="Envolventes", state="disabled")   

    def create_floorresp_action_frame(self):
        self.floorrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.floorrespFrame, text="Pisos", state="disabled")   
    
    def create_animation_action_frame(self):
        self.animationFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.animationFrame, text="Animacion", state="disabled")               
    
    def prepare_new_building(self):
        self.buildingName = tkSimpleDialog.askstring("Ingreso", "Ingresa el nombre del Edificio:", parent=self)
        self.prepare_to_building()
        #Disable sim and plot tabs
        self.disabled_tabs()

    def load_building(self):
        filename = askopenfilename(title="Elige un archivo para abrir", initialdir=".", filetypes=[("BuildingSim File","*.bsim")])
        tree = ET.parse(filename)
        root = tree.getroot()
        self.buildingName = root.get("name")        
        self.prepare_to_building()
        
        self.nameLabel.config(text=self.buildingName)
        
        self.dxEntry.insert(0, root.get("dx"))
        self.dyEntry.insert(0, root.get("dy"))
        
        for floorNode in tree.getiterator("floor"):
            floor = {"h":floorNode.get("h"), "E":floorNode.get("E"), "I":floorNode.get("I"), "gamma":floorNode.get("gamma"), "esp":floorNode.get("esp")}
            self.buildings.append(floor)
            self.update_floor_listbox()
            
        #Disable sim and plot tabs
        self.disabled_tabs()
            
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

    def create_styles(self):
        #Create styles         
        self.mainBg = "white"        
        style = Style()        
        style.configure("TLabel", font=("Verdana", 12))
        style.configure("Title.TLabel", font=("Verdana", 14))
        style.configure("Red.TLabel", font=("Verdana", 14), foreground="red")
    
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
        PanedWindow.__init__(self, master, orient=HORIZONTAL, width=master.winfo_screenwidth(), height=master.winfo_screenheight())
        #Create frames in paned window
        self.grid(padx=20, pady=20, sticky=N+S+E+W)
        self.grid_propagate(0)
        self.frameLeft = Frame(self, width=750)
        self.frameRight = Frame(self)
        self.add(self.frameLeft)
        self.add(self.frameRight)
        
        #Config master parameters
        self.master.config(bg="#003C6E")        
        self.master.title("Calculadora de Comportamiento de Edificios")        
        
        self.create_styles()
        self.create_menu()
        
        #Display ING logo
        self.canvas = Canvas(width = 77, height = 98)        
        self.gifLogo = PhotoImage(file = 'imgs/ing.gif')
        self.canvas.create_image(0, 0, image = self.gifLogo, anchor = NW)
        self.canvas.grid(row=0, column=0, padx=0, pady=0, sticky=NW)
        
        #Init variables
        self.buildingName = ""
        self.isCreated = False

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