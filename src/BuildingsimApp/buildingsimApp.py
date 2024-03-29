#!/usr/bin/python

from Tkinter import *
import tkSimpleDialog
import tkMessageBox
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
from senales import *

class App(PanedWindow):
    def unit_change(self, event):
        self.change_units()
        
    def change_units(self):
        distanceUnit = self.distanceUnitList[self.distanceUnitCombobox.current()]
        forceUnit = self.forceUnitList[self.forceUnitCombobox.current()]
        
        self.labelDx.config(text="Largo [" + distanceUnit + "]")
        self.labelDy.config(text="Ancho [" + distanceUnit + "]")
        self.labelH.config(text="Altura [" + distanceUnit + "]")
        self.labelE.config(text="Modulo Elasticidad [" + forceUnit + "/" + distanceUnit + "^2]")
        self.labelI.config(text="Momento Inercia [" + distanceUnit + "^4]")
        self.labelGamma.config(text="Peso Unitario Losa [" + forceUnit + "/" + distanceUnit + "^3]")
        self.labelEsp.config(text="Espesor [" + distanceUnit + "]")
        
    def calculate_g(self, distanceUnit):
        g = 9.806
        if distanceUnit == "cm":
            g = 980.6
        elif distanceUnit == "mm":
            g = 9806.
            
        return g
    
    def prepare_action(self):
        self.building = {}
        self.building['dx'] = float(self.dxEntry.get())
        self.building['dy'] =  float(self.dyEntry.get())
        self.building['g'] = self.calculate_g(self.distanceUnitList[self.distanceUnitCombobox.current()])
        self.building['xsi'] = 5. 				# Modal damping  [%]

        self.building['h'] = zeros((len(self.buildings),))
        self.building['E'] = zeros((len(self.buildings),))
        self.building['I'] = zeros((len(self.buildings),))
        self.building['gamma'] = zeros((len(self.buildings),))
        self.building['esp'] = zeros((len(self.buildings),))

        for i, item in enumerate(self.buildings):
            for key,value in item.iteritems():
                self.building[key][i] = float(value)
        
        self.building['name'] = self.buildingName
                
        self.building, isstable = form(self.building, self.geomEff.get())
        
        if isstable != 0:
            #Enable sim and plot tabs
            self.actionsTabs.tab(1, state="normal")
            self.actionsTabs.tab(2, state="normal")        
            self.actionsTabs.tab(3, state="normal")
            self.actionsTabs.tab(4, state="normal")
            
            #Fill parameters        
            self.plotmodeModeCombobox.configure(values=arange(1,len(self.buildings)+1,1).tolist())
        else:
            tkMessageBox.showerror("Error", "Tu Edificio colapso, no es posible realizar mas calculos con el")
    
    def plot_freqresp(self):
        variable = self.freqrespVariableCombobox.current()
        if variable == -1:
            variable = 0
        
        try:
            f0 = float(self.freqrespF0Entry.get())
        except:
            f0 = 0.
        try:
            f1 = float(self.freqrespF1Entry.get())
        except:
            f1 = 10.
        try:
            nfreq = float(self.freqrespNfreqEntry.get())
        except:
            nfreq = 200.
        if self.freqrespPlottypeCombobox.current() != -1:
            plottype = self.plottypeDictionary[self.freqrespPlottypeCombobox.get()]
        else:
            plottype = "plot"
        
        handle = freqresp(self.building, variable, f0, f1, nfreq, plottype)
        pl.show()
        
    def plot_plotmode(self):        
        print self.plotmodeModeCombobox.current()
        sys.stdout.flush()
        if self.plotmodeModeCombobox.current() >= 0:
            mode = int(self.plotmodeModeCombobox.get())            
            plotmode(self.building, mode)
            pl.show()
        else:
            tkMessageBox.showerror("Error", "Debes elegir un modo")    
        
    def compare_buildings(self):
        if len(self.compareFileEntry.get()) == 0:
            tkMessageBox.showerror("Error", "Debes seleccionar un archivo con otro edificio para comparar")
        else:
            #Load and prepare other building
            otherfilename = self.compareFileEntry.get()                        
            otherbuilding, dUnit, fUnit = self.load_building_to_dictionary(otherfilename)
            otherbuilding, isstable = form(otherbuilding)
            
            #Load parameters
            variable = self.compareVariableCombobox.current()
            if variable == -1:
                variable = 0
            
            try:
                f0 = float(self.compareF0Entry.get())
            except:
                f0 = 0.
            try:
                f1 = float(self.compareF1Entry.get())
            except:
                f1 = 10.
            try:
                nfreq = float(self.compareNfreqEntry.get())
            except:
                nfreq = 200.
            if self.comparePlottypeCombobox.current() != -1:
                plottype = self.plottypeDictionary[self.comparePlottypeCombobox.get()]
            else:
                plottype = "plot"
            
            handles = compare_buildings(self.building, otherbuilding, variable, f0, f1, nfreq, plottype)            
            pl.show()
        
    def load_building_to_compare(self):
        filename = askopenfilename(title="Elige un archivo para el otro Edificio", initialdir=".", filetypes=[("BuildingSim File","*.bsim")])
        self.compareFileEntry.insert(0, filename)
    
    def create_input(self):
        siminput = {}
        try:
            tmax = float(self.simTmaxEntry.get())
        except:
            tmax = 60
        try:
            dt = float(self.simDtEntry.get())
        except:
            dt = 0.01
        
        t = arange(0., tmax , dt)
        siminput['t'] = t
        
        if self.simInputTypeCombobox.current() == -1:
            inputType = "armonica"
        else:
            inputType = self.simInputTypeCombobox.get()
            
        if inputType == "armonica":
            try:
                a = float(self.inputHarmonicAEntry.get())
            except:
                a = 9.806
            try:
                f = float(self.inputHarmonicFEntry.get())
            except:
                f = 1
            try:
                b = float(self.inputHarmonicBEntry.get())
            except:
                b = 0
            try:
                K = float(self.inputHarmonicKEntry.get())
            except:
                K = 0.
            try:
                tmin = float(self.inputHarmonicTminEntry.get())
            except:
                tmin = 0.
            try:
                tmax = float(self.inputHarmonicTmaxEntry.get())
            except:
                tmax = Inf
            siminput['ug'] = harmonic(t, a, f, b, K, tmin, tmax)
        elif inputType == "rectangular":
            try:
                a = float(self.inputRectangularAEntry.get())
            except:
                a = 9.806
            try:
                t0 = float(self.inputRectangularT0Entry.get())
            except:
                t0 = 1
            try:
                f = float(self.inputRectangularFEntry.get())
            except:
                f = 1
            try:
                tmin = float(self.inputRectangularTminEntry.get())
            except:
                tmin = 0.
            try:
                tmax = float(self.inputRectangularTmaxEntry.get())
            except:
                tmax = Inf
            siminput['ug'] = rectangular(t, a, t0, f, tmin, tmax)
        else: # TODO: Carga desde archivo .CSV
            siminput['ug'] = rectangular(t, 9.806)
        
        return siminput
        
    def show_input(self):
        siminput = self.create_input()
        plotinput(siminput)
        pl.show()
        
    def simulate(self):
        self.siminput = self.create_input()
        self.buildingSim = response(self.building, self.siminput)        
        
        #Enable animate and other plot tabs
        self.actionsTabs.tab(5, state="normal")
        self.actionsTabs.tab(6, state="normal")        
        self.actionsTabs.tab(7, state="normal")
        #Fill parameters        
        self.floorrespFloorCombobox.configure(values=arange(1,len(self.buildings)+1,1).tolist())
        
    def plot_envresp(self):
        if self.envrespVariableCombobox.current() == -1:
            variable = "acc"
        else:
            variable = self.envrespvariableDictionary[self.envrespVariableCombobox.get()]            
        envresp(self.building, self.siminput, self.buildingSim, variable)
        pl.show()
        
    def plot_floorresp(self):
        if self.floorrespFloorCombobox.current() >= 0:
            floornumber = int(self.floorrespFloorCombobox.get())
            floorresp(self.building, self.siminput, self.buildingSim, floornumber)
            pl.show()
        else:
            tkMessageBox.showerror("Error", "Debes elegir un piso")                    
        
    def animate_simulation(self):
        animdef(self.building, self.buildingSim['dis'])
    
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
                s = s + "[" + key + "=" + str(value) + "]"
            self.floorListBox.insert(END, s)
    
    def sim_input_type_change(self, event):
        inputType = self.simInputTypeCombobox.get()
        self.inputHarmonicFrame.grid_forget()
        self.inputRectangularFrame.grid_forget()
        if inputType == "armonica":
            self.inputHarmonicFrame.grid(row=4, column=0, columnspan=2, padx=3, pady=5, sticky=W)
        elif inputType == "rectangular":
            self.inputRectangularFrame.grid(row=4, column=0, columnspan=2, padx=3, pady=5, sticky=W)
        else: # TODO: Carga desde archivo .CSV
            print "other"
    
    def prepare_to_building(self):
        self.buildings = []
        if len(self.buildingName) != 0 and self.isCreated == False:
            self.isCreated = True            
            self.create_building_frame()
            self.create_actions_frame()
        else:
            self.dxEntry.delete(0, END)
            self.dyEntry.delete(0, END)
            self.floorListBox.delete(0, END)
    
    def create_building_frame(self):
        currentRow = 0
        #Create Widgets            
        Label(self.frameLeft, text="Edificio ", style="Title.TLabel").grid(row=currentRow, column=0, padx=5, pady=5, sticky=NE)
        self.nameLabel = Label(self.frameLeft, text=self.buildingName, style="Red.TLabel")
        self.nameLabel.grid(row=currentRow, column=1, padx=5, pady=5, sticky=NW)
        
        # Unidades
        currentRow += 1
        Label(self.frameLeft, text="Unidades", style="Title.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=10)
        currentRow += 1
        Label(self.frameLeft, text="Al cambiar las unidades los valores NO se cambian automaticamente", style="VerySmall.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=10)
        currentRow += 1
        Label(self.frameLeft, text="Distancia").grid(row=currentRow, column=0, padx=5, sticky=E)        
        self.distanceUnitList = ("m", "cm", "mm")
        self.distanceUnitCombobox = Combobox(self.frameLeft, values=("metros [m]", "centimetros [cm]", "milimetros [mm]"))
        self.distanceUnitCombobox.grid(row=currentRow, column=1, padx=3, pady=5, sticky=W)
        self.distanceUnitCombobox.current(newindex=0)
        self.distanceUnitCombobox.bind('<<ComboboxSelected>>', self.unit_change)
        currentRow += 1
        Label(self.frameLeft, text="Fuerza").grid(row=currentRow, column=0, padx=5, sticky=E)        
        self.forceUnitList = ("tonf", "Kgf", "N", "KN")
        self.forceUnitCombobox = Combobox(self.frameLeft, values=("toneladas fuerza [tonf]", "kilos fuerza [Kgf]", "Newton [N]", "Kilonewton [KN]"))
        self.forceUnitCombobox.grid(row=currentRow, column=1, padx=3, pady=5, sticky=W)
        self.forceUnitCombobox.current(newindex=0)
        self.forceUnitCombobox.bind('<<ComboboxSelected>>', self.unit_change)

        # Base
        currentRow += 1
        Label(self.frameLeft, text="Base", style="Title.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=10)
        currentRow += 1
        self.labelDx = Label(self.frameLeft, text="Largo [m]")
        self.labelDx.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.dxEntry = Entry(self.frameLeft)
        self.dxEntry.grid(row=currentRow, column=1, padx=3)
        currentRow += 1
        self.labelDy = Label(self.frameLeft, text="Ancho [m]")
        self.labelDy.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.dyEntry = Entry(self.frameLeft)
        self.dyEntry.grid(row=currentRow, column=1, padx=3)
        
        currentRow += 1

        # Pisos
        Label(self.frameLeft, text="Pisos", style="Title.TLabel").grid(row=currentRow, column=0, columnspan=3, padx=5, pady=5)
        currentRow += 1
        self.labelH = Label(self.frameLeft, text="Altura [m]")
        self.labelH.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.hEntry = Entry(self.frameLeft)
        self.hEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        self.labelE = Label(self.frameLeft, text="Modulo Elasticidad [tonf/m^2]")
        self.labelE.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.eEntry = Entry(self.frameLeft)
        self.eEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        self.labelI = Label(self.frameLeft, text="Momento Inercia [m^4]")
        self.labelI.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.iEntry = Entry(self.frameLeft)
        self.iEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        self.labelGamma = Label(self.frameLeft, text="Peso Unitario Losa [tonf/m^3]")
        self.labelGamma.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.gammaEntry = Entry(self.frameLeft)
        self.gammaEntry.grid(row=currentRow, column=1, padx=5)
        
        currentRow += 1
        self.labelEsp = Label(self.frameLeft, text="Espesor [m]")
        self.labelEsp.grid(row=currentRow, column=0, padx=5, sticky=E)
        self.espEntry = Entry(self.frameLeft)
        self.espEntry.grid(row=currentRow, column=1, padx=5)
                    
        self.addFloorButton = Button(self.frameLeft, text="Agregar >>", command=self.add_floor, width=10)
        self.addFloorButton.grid(row=10, column=2, rowspan=2, padx=5)

        self.changeFloorButton = Button(self.frameLeft, text="Cambiar >>", command=self.change_floor, width=10)
        self.changeFloorButton.grid(row=11, column=2, rowspan=2, padx=5)
        
        self.editFloorButton = Button(self.frameLeft, text="<< Editar", command=self.edit_floor, width=10)
        self.editFloorButton.grid(row=12, column=2, rowspan=2, padx=5)

        self.floorListBox = Listbox(self.frameLeft, width=40)
        self.floorListBox.grid(row=10, column=3, rowspan=4, columnspan=2, padx=10)

        self.deleteFloorButton = Button(self.frameLeft, text="Eliminar", command=self.delete_floor, width=10)
        self.deleteFloorButton.grid(row=14, column=4, padx=5)

    def create_actions_frame(self):
        self.actionsTabs = Notebook(self.frameRight)
        self.actionsTabs.pack(fill=BOTH, expand=1)
        self.create_prepare_action_frame()
        self.create_freqresp_action_frame()
        self.create_plotmode_action_frame()                
        self.create_compare_action_frame()
        self.create_sim_action_frame()
        self.create_envresp_action_frame()
        self.create_floorresp_action_frame()
        self.create_animation_action_frame()
    
    def create_prepare_action_frame(self):
        self.prepareFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.prepareFrame, text="Preparar")
        
        Label(self.prepareFrame, text="Preparar Calculos Edificio", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        self.geomEff = IntVar()
        self.geomEffCheckbox = Checkbutton(self.prepareFrame, text="Incluir efectos de 2do orden", variable=self.geomEff)
        self.geomEffCheckbox.grid(row=1, column=0, columnspan=2, padx=3, pady=5, sticky=W)
        
        # Prepare Button
        self.prepareButton = Button(self.prepareFrame, text="Preparar", command=self.prepare_action)
        self.prepareButton.grid(row=3, column=0, padx=10, pady=10)
        
    def create_sim_action_frame(self):
        self.simFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.simFrame, text="Simular", state="disabled")
        
        Label(self.simFrame, text="Simular", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        Label(self.simFrame, text="Tiempo de Simulacion [s]", style="Small.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.simTmaxEntry = Entry(self.simFrame)
        self.simTmaxEntry.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.simFrame, text="Intervalo [s]", style="Small.TLabel").grid(row=2, column=0, padx=5, sticky=E)        
        self.simDtEntry = Entry(self.simFrame)
        self.simDtEntry.grid(row=2, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.simFrame, text="Tipo de Movimiento", style="Small.TLabel").grid(row=3, column=0, padx=5, sticky=E)        
        self.simInputTypeCombobox = Combobox(self.simFrame, values=("armonica", "rectangular"))
        self.simInputTypeCombobox.grid(row=3, column=1, padx=3, pady=5, sticky=W)
        self.simInputTypeCombobox.bind('<<ComboboxSelected>>', self.sim_input_type_change)
        
        # View input Button
        self.simShowInputButton = Button(self.simFrame, text="Ver Movimiento", command=self.show_input)
        self.simShowInputButton.grid(row=9, column=1, padx=10, pady=10)
        
        # Plot Button
        self.simButton = Button(self.simFrame, text="Simular", command=self.simulate)
        self.simButton.grid(row=10, column=0, padx=10, pady=10)
        
        # Create frames for different inputs
        # Harmonic
        self.inputHarmonicFrame = Frame(self.simFrame)
        Label(self.inputHarmonicFrame, text="Parametros Armonico", style="Small.TLabel").grid(row=0, column=0, padx=5, sticky=W+E+N+S)
        
        Label(self.inputHarmonicFrame, text="Amplitud A [m/s^2]", style="VerySmall.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.inputHarmonicAEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicAEntry.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputHarmonicFrame, text="Frecuencia [Hz]", style="VerySmall.TLabel").grid(row=2, column=0, padx=5, sticky=E)        
        self.inputHarmonicFEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicFEntry.grid(row=2, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputHarmonicFrame, text="Amplitud B [m/s^2]", style="VerySmall.TLabel").grid(row=3, column=0, padx=5, sticky=E)        
        self.inputHarmonicBEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicBEntry.grid(row=3, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputHarmonicFrame, text="Tasa Crecimiento Frecuencia [Hz/s]", style="VerySmall.TLabel").grid(row=4, column=0, padx=5, sticky=E)        
        self.inputHarmonicKEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicKEntry.grid(row=4, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputHarmonicFrame, text="Tiempo Inicio [s]", style="VerySmall.TLabel").grid(row=5, column=0, padx=5, sticky=E)        
        self.inputHarmonicTminEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicTminEntry.grid(row=5, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputHarmonicFrame, text="Tiempo Final [s]", style="VerySmall.TLabel").grid(row=6, column=0, padx=5, sticky=E)        
        self.inputHarmonicTmaxEntry = Entry(self.inputHarmonicFrame)
        self.inputHarmonicTmaxEntry.grid(row=6, column=1, padx=3, pady=5, sticky=W)
        
        # Rectangular
        self.inputRectangularFrame = Frame(self.simFrame)
        Label(self.inputRectangularFrame, text="Parametros Rectangular", style="Small.TLabel").grid(row=0, column=0, padx=5, sticky=W+E+N+S)
        
        Label(self.inputRectangularFrame, text="Amplitud Aceleracion [m/s^2]", style="VerySmall.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.inputRectangularAEntry = Entry(self.inputRectangularFrame)
        self.inputRectangularAEntry.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputRectangularFrame, text="Frecuencia [Hz]", style="VerySmall.TLabel").grid(row=2, column=0, padx=5, sticky=E)        
        self.inputRectangularFEntry = Entry(self.inputRectangularFrame)
        self.inputRectangularFEntry.grid(row=2, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputRectangularFrame, text="Desfase [s]", style="VerySmall.TLabel").grid(row=3, column=0, padx=5, sticky=E)        
        self.inputRectangularT0Entry = Entry(self.inputRectangularFrame)
        self.inputRectangularT0Entry.grid(row=3, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputRectangularFrame, text="Tiempo Inicio [s]", style="VerySmall.TLabel").grid(row=4, column=0, padx=5, sticky=E)        
        self.inputRectangularTminEntry = Entry(self.inputRectangularFrame)
        self.inputRectangularTminEntry.grid(row=4, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.inputRectangularFrame, text="Tiempo Final [s]", style="VerySmall.TLabel").grid(row=5, column=0, padx=5, sticky=E)        
        self.inputRectangularTmaxEntry = Entry(self.inputRectangularFrame)
        self.inputRectangularTmaxEntry.grid(row=5, column=1, padx=3, pady=5, sticky=W)
        
    def create_freqresp_action_frame(self):
        self.freqrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.freqrespFrame, text="Resp. a Frecuencia", state="disabled")
        
        Label(self.freqrespFrame, text="Grafico de Respuesta a Frecuencia", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        Label(self.freqrespFrame, text="Variable de Salida", style="Small.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.freqrespVariableCombobox = Combobox(self.freqrespFrame, values=("desplazamiento", "velocidad", "aceleracion"))
        self.freqrespVariableCombobox.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.freqrespFrame, text="Frecuencia Minima", style="Small.TLabel").grid(row=2, column=0, padx=5, sticky=E)        
        self.freqrespF0Entry = Entry(self.freqrespFrame)
        self.freqrespF0Entry.grid(row=2, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.freqrespFrame, text="Frecuencia Maxima", style="Small.TLabel").grid(row=3, column=0, padx=5, sticky=E)        
        self.freqrespF1Entry = Entry(self.freqrespFrame)
        self.freqrespF1Entry.grid(row=3, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.freqrespFrame, text="Numero de Frecuencias", style="Small.TLabel").grid(row=4, column=0, padx=5, sticky=E)        
        self.freqrespNfreqEntry = Entry(self.freqrespFrame)
        self.freqrespNfreqEntry.grid(row=4, column=1, padx=3, pady=5, sticky=W)
        
        self.plottypeDictionary = {"simple":"plot", "logaritmico doble":"loglog", "semilogaritmico en x":"semilogx", "semilogaritmico en y":"semilogy"}
        Label(self.freqrespFrame, text="Tipo de Funcion", style="Small.TLabel").grid(row=5, column=0, padx=5, sticky=E)        
        self.freqrespPlottypeCombobox = Combobox(self.freqrespFrame, values=("simple", "logaritmico doble", "semilogaritmico en x", "semilogaritmico en y"))
        self.freqrespPlottypeCombobox.grid(row=5, column=1, padx=3, pady=5, sticky=W)
        
        # Plot Button
        self.freqrespButton = Button(self.freqrespFrame, text="Graficar", command=self.plot_freqresp)
        self.freqrespButton.grid(row=8, column=0, padx=10, pady=10)
        
    def create_plotmode_action_frame(self):
        self.plotmodeFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.plotmodeFrame, text="Modos", state="disabled")   
        
        Label(self.plotmodeFrame, text="Mostrar modos y periodos naturales de vibracion.", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        Label(self.plotmodeFrame, text="Modo", style="Small.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.plotmodeModeCombobox = Combobox(self.plotmodeFrame)
        self.plotmodeModeCombobox.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        # Plot Button
        self.plotmodeButton = Button(self.plotmodeFrame, text="Graficar", command=self.plot_plotmode)
        self.plotmodeButton.grid(row=8, column=0, padx=10, pady=10)
        
    def create_compare_action_frame(self):
        self.compareFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.compareFrame, text="Comparar", state="disabled")
        
        Label(self.compareFrame, text="Comparar Edificaciones", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        Label(self.compareFrame, text="Otro Edificio").grid(row=1, column=0, columnspan=2, padx=5, sticky=N+S+W+E)        
        self.compareFileEntry = Entry(self.compareFrame)
        self.compareFileEntry.grid(row=2, column=0, padx=3, pady=5, sticky=W)        
        self.compareLoadFileButton = Button(self.compareFrame, text="Seleccionar Archivo", command=self.load_building_to_compare)
        self.compareLoadFileButton.grid(row=2, column=1, padx=10, pady=10)
        
        Label(self.compareFrame, text="Parametros para Comparacion").grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)        
        Label(self.compareFrame, text="Variable de Salida", style="Small.TLabel").grid(row=4, column=0, padx=5, sticky=E)        
        self.compareVariableCombobox = Combobox(self.compareFrame, values=("desplazamiento", "velocidad", "aceleracion"))
        self.compareVariableCombobox.grid(row=4, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.compareFrame, text="Frecuencia Minima", style="Small.TLabel").grid(row=5, column=0, padx=5, sticky=E)        
        self.compareF0Entry = Entry(self.compareFrame)
        self.compareF0Entry.grid(row=5, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.compareFrame, text="Frecuencia Maxima", style="Small.TLabel").grid(row=6, column=0, padx=5, sticky=E)        
        self.compareF1Entry = Entry(self.compareFrame)
        self.compareF1Entry.grid(row=6, column=1, padx=3, pady=5, sticky=W)
        
        Label(self.compareFrame, text="Numero de Frecuencias", style="Small.TLabel").grid(row=7, column=0, padx=5, sticky=E)        
        self.compareNfreqEntry = Entry(self.compareFrame)
        self.compareNfreqEntry.grid(row=7, column=1, padx=3, pady=5, sticky=W)
                
        Label(self.compareFrame, text="Tipo de Funcion", style="Small.TLabel").grid(row=8, column=0, padx=5, sticky=E)        
        self.comparePlottypeCombobox = Combobox(self.compareFrame, values=("simple", "logaritmico doble", "semilogaritmico en x", "semilogaritmico en y"))
        self.comparePlottypeCombobox.grid(row=8, column=1, padx=3, pady=5, sticky=W)
        
        # Compare Button
        self.compareButton = Button(self.compareFrame, text="Comparar", command=self.compare_buildings)
        self.compareButton.grid(row=12, column=0, padx=10, pady=10)
    
    def create_envresp_action_frame(self):
        self.envrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.envrespFrame, text="Envolventes", state="disabled")   
        
        Label(self.envrespFrame, text="Grafico de Envolventes", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        self.envrespvariableDictionary = {"aceleracion":"acc", "velocidad":"vel", "desplazamiento":"dis", "deriva":"dri"}
        Label(self.envrespFrame, text="Variable de Salida", style="Small.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.envrespVariableCombobox = Combobox(self.envrespFrame, values=("aceleracion", "velocidad", "desplazamiento", "deriva"))
        self.envrespVariableCombobox.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        # Plot Button
        self.envrespButton = Button(self.envrespFrame, text="Graficar", command=self.plot_envresp)
        self.envrespButton.grid(row=5, column=0, padx=10, pady=10)

    def create_floorresp_action_frame(self):
        self.floorrespFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.floorrespFrame, text="Pisos", state="disabled")   
        
        Label(self.floorrespFrame, text="Grafico Respuestas por Piso", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        Label(self.floorrespFrame, text="Piso", style="Small.TLabel").grid(row=1, column=0, padx=5, sticky=E)        
        self.floorrespFloorCombobox = Combobox(self.floorrespFrame)
        self.floorrespFloorCombobox.grid(row=1, column=1, padx=3, pady=5, sticky=W)
        
        # Plot Button
        self.floorrespButton = Button(self.floorrespFrame, text="Graficar", command=self.plot_floorresp)
        self.floorrespButton.grid(row=8, column=0, padx=10, pady=10)
    
    def create_animation_action_frame(self):
        self.animationFrame = Frame(self.actionsTabs)
        self.actionsTabs.add(self.animationFrame, text="Animacion", state="disabled")               
        
        Label(self.animationFrame, text="Visualizar Animacion de Simulacion", style="Title.TLabel").grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky=N+S+W+E)
        
        # Animation Button
        self.animationButton = Button(self.animationFrame, text="Animar", command=self.animate_simulation)
        self.animationButton.grid(row=1, column=0, padx=10, pady=10)
    
    def prepare_new_building(self):
        self.buildingName = tkSimpleDialog.askstring("Ingreso", "Ingresa el nombre del Edificio:", parent=self)
        self.prepare_to_building()
        #Disable sim and plot tabs
        self.disabled_tabs()
    
    def load_building_to_dictionary(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        
        building = {}
        building['name'] = root.get("name")
        building['dx'] = float(root.get("dx"))
        building['dy'] =  float(root.get("dy"))        
        
        distanceUnit = root.get("du", "m")
        forceUnit = root.get("fu", "tonf")
        
        building['du'] = distanceUnit
        building['fu'] = forceUnit
        building['g'] = self.calculate_g(distanceUnit)
        building['xsi'] = 5. 				# Modal damping  [%]        
                
        numfloors = sum(1 for floor in tree.getiterator("floor"))
        
        building['h'] = zeros((numfloors,))       # Altura entrepiso [m]
        building['E'] = zeros((numfloors,))       # Modulo elasticidad piso [tonf/cm**2]
        building['I'] = zeros((numfloors,))       # Momento inercia columnas [m**4]
        building['gamma'] = zeros((numfloors,))   # Peso unitario losas [tonf/m**3]
        building['esp'] = zeros((numfloors,))     # espesor losas [m]
        
        for i, floorNode in enumerate(tree.getiterator("floor")):
            floor = {"h":floorNode.get("h"), "E":floorNode.get("E"), "I":floorNode.get("I"), "gamma":floorNode.get("gamma"), "esp":floorNode.get("esp")}
            for key,value in floor.iteritems():
                building[key][i] = float(value) 
        
        return building, distanceUnit, forceUnit
        
    def load_building(self):
        filename = askopenfilename(title="Elige un archivo para abrir", initialdir=".", filetypes=[("BuildingSim File","*.bsim")])
        building, dUnit, fUnit = self.load_building_to_dictionary(filename)
                
        self.buildingName = building['name']
        self.prepare_to_building()
        
        self.distanceUnitCombobox.current(newindex=self.distanceUnitList.index(dUnit))
        self.forceUnitCombobox.current(newindex=self.forceUnitList.index(fUnit))
        self.change_units()
        
        self.nameLabel.config(text=self.buildingName)
        
        self.dxEntry.insert(0, building['dx'])
        self.dyEntry.insert(0, building['dy'])
        
        self.buildings = []
        for i, value in enumerate(building['h']):            
            floor = {"h":building["h"][i], "E":building["E"][i], "I":building["I"][i], "gamma":building["gamma"][i], "esp":building["esp"][i]}
            self.buildings.append(floor)
            
        self.update_floor_listbox()
            
        #Disable sim and plot tabs
        self.disabled_tabs()
            
    def save_building(self):
        filename = asksaveasfilename(parent=self.master, title="Ingresa o selecciona el archivo donde guardar", filetypes=[("Buildingsim File","*.bsim")])

        buildingNode = ET.Element("building")
        buildingNode.set("name", self.buildingName)
        
        distanceUnit = self.distanceUnitList[self.distanceUnitCombobox.current()]
        forceUnit = self.forceUnitList[self.forceUnitCombobox.current()]        
        buildingNode.set("du", distanceUnit)
        buildingNode.set("fu", forceUnit)
        
        buildingNode.set("dx", self.dxEntry.get())
        buildingNode.set("dy", self.dyEntry.get())

        floorsListNode = ET.SubElement(buildingNode, "floors")

        #Save floors
        for i, item in enumerate(self.buildings):
            floorNode = ET.SubElement(floorsListNode, "floor")
            for key,value in item.iteritems():
                floorNode.set(key, str(value))

        tree = ET.ElementTree(buildingNode)
        tree.write(filename)

    def create_styles(self):
        #Create styles         
        self.mainBg = "white"        
        style = Style()        
        style.configure("TLabel", font=("Verdana", 12))
        style.configure("Title.TLabel", font=("Verdana", 14))
        style.configure("Red.TLabel", font=("Verdana", 14), foreground="red")
        style.configure("Small.TLabel", font=("Verdana", 10))
        style.configure("VerySmall.TLabel", font=("Verdana", 9))
    
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
