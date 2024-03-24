import customtkinter as ctk
import math
import numpy as np
import random

class Torus(ctk.CTk):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.canvassize = [600,600]
        self.maincanvas = ctk.CTkCanvas(self, width = self.canvassize[0], height = self.canvassize[1], bg = '#181818')
        self.maincanvas.pack(fill = ctk.BOTH, expand=True)
        self.colours = self.getcolours()
        self.bindings()
        self.initialize

    def initialize(self):
        self.maincanvas.delete("all")
        self.origin = self.maincanvas.create_aa_circle(math.floor(self.canvassize[0]/2),math.floor(self.canvassize[1]/2),fill="white",radius = 1)
        self.tris = []
        self.points = self.parameterize()
        self.drawsurface()

    def bindings(self):
        self.maincanvas.bind("<Motion>", self.getcursorpos)
        self.maincanvas.bind("<ButtonPress-1>", lambda event: self.userrotates(np.array((event.x-self.canvassize[0]/2,event.y-self.canvassize[1]/2)), True, event.widget.find_withtag(ctk.CURRENT)))
        self.maincanvas.bind("<ButtonPress-3>", lambda event: self.usertranslates(np.array((event.x-self.canvassize[0]/2,event.y-self.canvassize[1]/2)), True))
        self.maincanvas.bind("<ButtonRelease-1>", self.button1release)
        self.maincanvas.bind("<ButtonRelease-3>", self.button3release)
        self.maincanvas.bind("<MouseWheel>", self.on_scroll)
        self.maincanvas.bind("<Configure>", self.on_canvas_resize)
        self.maincanvas.bind("<Double-Button-1>", lambda event: self.initialize())

    def on_canvas_resize(self,event):
        self.canvassize = [event.width,event.height]
        self.initialize()

    def on_scroll(self, event):
        if event.delta < 0:
            scale = 0.95
        if event.delta > 0:
            scale = 1/0.95

        self.points = [[coord*scale for coord in point] for point in self.points]
        self.rotatesurface(0, 0, 0)

    def button1release(self, event):
        self.button1state = 0

    def button3release(self,event):
        self.button3state = 0

    def getcursorpos(self, event):
        self.cursorpos = np.array((event.x-self.canvassize[0]/2,event.y-self.canvassize[1]/2))

    def usertranslates(self, cursorposatcall, firstcall):
        try:
            self.after_cancel(self.tranlationloop)
        except AttributeError:
            pass
        if firstcall:
            self.button3state = 1
        displacement = self.cursorpos - cursorposatcall
        displacementwithz = np.array((displacement[0],displacement[1],0))
        self.points = [point+displacementwithz for point in self.points]
        self.rotatesurface(0,0,0)
        if self.button3state == 1:
            self.after(50, self.usertranslates, self.cursorpos, False)
        elif np.linalg.norm(displacement)!=0:
            self.translatewithmomentum(displacementwithz)

    def translatewithmomentum(self, displacement):
        self.points = [point+displacement for point in self.points]
        self.rotatesurface(0,0,0)
        self.tranlationloop = self.after(50, self.translatewithmomentum, displacement*0.8)

    def userrotates(self, cursorposatcall, firstcall, item):
        try:
            self.after_cancel(self.rotationloop)
        except AttributeError:
            pass
        if firstcall:
            self.button1state = 1
        displacement = cursorposatcall-self.cursorpos
        if np.linalg.norm(displacement) != 0:
            cursordirection = self.cursorpos/np.linalg.norm(self.cursorpos)
            xyrotation = abs(np.dot(displacement,cursordirection))/np.linalg.norm(displacement)
            thetaz = -np.cross(self.cursorpos,displacement)/(np.linalg.norm(displacement))/3000
            thetax = displacement[1]*math.pi/self.canvassize[1]*xyrotation
            thetay = -displacement[0]*math.pi/self.canvassize[0]*xyrotation
            self.rotatesurface(thetax, thetay, thetaz)

        if self.button1state == 1:
            self.after(50, self.userrotates, self.cursorpos, False, item)
        elif np.linalg.norm(displacement) != 0:
            self.rotatewithmomentum(thetax,thetay,thetaz)
    
    def rotatewithmomentum(self, thetax, thetay, thetaz):
        self.rotatesurface(thetax,thetay,thetaz)
        self.rotationloop = self.after(50, self.rotatewithmomentum,thetax,thetay,thetaz)

    def rotatesurface(self, thetax, thetay, thetaz):
        trisbyz = sorted(self.tri_indicies, key = self.sorter)
        for i, point in enumerate(self.points):
            self.points[i] = self.Rx(thetax) @ self.Ry(thetay) @ self.Rz(thetaz) @ np.transpose(point)       # type: ignore
        for tri_index in trisbyz:
            i = self.tri_indicies.index(tri_index)
            self.maincanvas.coords(self.tris[i], self.points[tri_index[0]][0]+0.5*self.canvassize[0],self.points[tri_index[0]][1]+0.5*self.canvassize[1],
                                    self.points[tri_index[1]][0]+0.5*self.canvassize[0],self.points[tri_index[1]][1]+0.5*self.canvassize[1],
                                    self.points[tri_index[2]][0]+0.5*self.canvassize[0],self.points[tri_index[2]][1]+0.5*self.canvassize[1])
            self.maincanvas.lift(self.tris[i])
            
            averagez = 1/3*(self.points[tri_index[0]][2]+self.points[tri_index[1]][2]+self.points[tri_index[2]][2])
            if averagez<0:
                self.maincanvas.lift(self.origin)
            colourindex = math.floor(np.interp(averagez,np.linspace(-0.3*self.canvassize[0],0.3*self.canvassize[0],len(self.colours)),np.arange(len(self.colours))[::-1]))
            self.maincanvas.itemconfig(self.tris[i], fill = self.colours[colourindex])

    def sorter(self, tri_index):
        isshown = self.maincanvas.itemcget(self.tris[self.tri_indicies.index(tri_index)],"state") == "normal"
        averagez = 1/3*(self.points[tri_index[0]][2]+self.points[tri_index[1]][2]+self.points[tri_index[2]][2])
        return (not isshown, averagez)

    def drawsurface(self):
        delay = 0
        for tri_index in self.tri_indicies:
            self.drawtriangle(tri_index)
        self.rotatewithmomentum(0.02,0.02,0.02)
        for tri in self.tris:
            self.after(delay,self.fillintri,tri)
            delay += 3  
    
    def fillintri(self, tri):
        self.maincanvas.itemconfig(tri, state="normal")
        self.maincanvas.lower(tri)

    def drawtriangle(self, tri_index):
        averagez = 1/3*(self.points[tri_index[0]][2]+self.points[tri_index[1]][2]+self.points[tri_index[2]][2])
        
        colourindex = math.floor(np.interp(averagez,np.linspace(-0.3*self.canvassize[0],0.3*self.canvassize[0],len(self.colours)),np.arange(len(self.colours))[::-1]))
        triangle = self.maincanvas.create_polygon(self.points[tri_index[0]][0]+0.5*self.canvassize[0],self.points[tri_index[0]][1]+0.5*self.canvassize[1],
                                                self.points[tri_index[1]][0]+0.5*self.canvassize[0],self.points[tri_index[1]][1]+0.5*self.canvassize[1],
                                                self.points[tri_index[2]][0]+0.5*self.canvassize[0],self.points[tri_index[2]][1]+0.5*self.canvassize[1],
                                                fill = self.colours[colourindex], state="hidden")
        self.tris.append(triangle)

    def getcolours(self):
        with open('blue gradient.txt','r') as file:
            gradient = file.read()
            gradient = gradient.replace('\n',' ')
            colours = gradient.split(' ')
            return colours

    def Rx(self, psi):
        return np.array([[1,0,0],
                         [0,math.cos(psi),-math.sin(psi)],
                         [0,math.sin(psi),math.cos(psi)]])
    
    def Ry(self, psi):
        return np.array([[math.cos(psi),0,math.sin(psi)],
                         [0,1,0],
                         [-math.sin(psi),0,math.cos(psi)]])
    
    def Rz(self, psi):
        return np.array([[math.cos(psi),-math.sin(psi),0],
                         [math.sin(psi),math.cos(psi),0],
                         [0,0,1]])

    def maketrigrid(self,ucount,vcount):
        trimap = []
        for i in range(vcount-1):
            for j in range(ucount-1):
                trimap.append([i+j*vcount,i+1+j*vcount,i+vcount+j*vcount])
                trimap.append([i+1+j*vcount,i+vcount+j*vcount,i+vcount+1+j*vcount])
        return trimap
    
    def torusparam(self,u,v,r,R):
        i = math.cos(u)*(R+math.cos(v)*r)
        j = math.sin(u)*(R+math.cos(v)*r)
        k = math.sin(v)
        return np.array((i, j, k))

    def mobiusstripparam(self,u,v):
        i = (1-v*math.cos(u/2))*math.cos(u)
        j = (1-v*math.cos(u/2))*math.sin(u)
        k = v*math.sin(u/2)
        return np.array((i, j, k))

    def surface(self, *args):
        return self.torusparam(*args)

    def getuvratio(self, max, ubounds, vbounds, *args):
        dudv = 0
        n = 5
        
        for i in range(n):
            un1 = random.uniform(ubounds[0], ubounds[1])
            vn1 = random.uniform(vbounds[0], vbounds[1])
            un2 = un1 + np.sign(un1-(sum(ubounds)/2-un1))*1/10*(ubounds[1]-ubounds[0])
            vn2 = vn1 + np.sign(un1-(sum(vbounds)/2-vn1))*1/10*(vbounds[1]-vbounds[0])
            du = np.linalg.norm(self.surface(un1, vn1, *args)-self.surface(un2,vn1, *args))
            dv = np.linalg.norm(self.surface(un1, vn1, *args)-self.surface(un1,vn2, *args))
            if du/dv:
                dudv += du/dv
        dudv *= 1/n
        vcount = math.floor(math.sqrt(1/dudv*max))
        ucount = math.ceil(dudv*vcount)
        return ucount, vcount

    def parameterize(self):
        littler = 1
        bigr = 3
        ubounds = [-math.pi,math.pi]
        vbounds = [0,2*math.pi]

        ucount, vcount = self.getuvratio(500,ubounds,vbounds,littler,bigr)
        min = float('inf')
        max = float('-inf')
        absolutepoints = []
        for u in np.linspace(ubounds[0],ubounds[1],ucount):
            for v in np.linspace(vbounds[0],vbounds[1],vcount):
                point = self.surface(u,v,littler,bigr)
                for i in range(3):
                    if point[i] < min:
                        min = point[i]
                    if point[i] > max:
                        max = point[i]
                absolutepoints.append(point)
        
        self.tri_indicies = self.maketrigrid(ucount,vcount)
        
        canvascoords = [[self.canvassize[1]*0.6/(max-min)*(point[i]-min)-self.canvassize[1]*0.3 for i in range(3)] for point in absolutepoints]
        return(canvascoords)

mainwindow = Torus()
mainwindow.mainloop()
