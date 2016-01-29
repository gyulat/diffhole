#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Diffraction through pinholes

Inspired by:
http://math.stackexchange.com/questions/733754/visually-stunning-math-concepts-which-are-easy-to-explain
(see no.36):
"Fourier transform of the light intensity due to a diffraction pattern caused by light going through 
 8 pinholes and interfering on a wall"
also https://www.youtube.com/watch?v=1UVbUWuyNmk by Math Doobler

Initial version:

- preview effect of various settings:
    - wavenumber
    - pattern shift (x and y)
    - number of pinholes
    - editable color palette
- saves images

Dependencies:

- PyQt4 (sourceforge.net/projects/pyqt/files/PyQt4/)
- pyqtgraph (www.pyqtgraph.org)
- numpy (www.numpy.org)
- scipy (www.scipy.org)

author: Gyula Tóth
last edited: January 2016
"""

import sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
from scipy import fftpack

## math of the image
def pindr(k=350, u=1.0, v=1.0, npx = 600, nholes=8, R=2.0, d=0.1):
    """
    FFT of diffraction pattern through 'npins' pinholes
    k      : wavenumber
    u, v   : shift of pattern (0 - 2pi) in x and y directions
    npx    : size of image in pixels
    nholes : number of pinholes
    R      : radius of pinholes' circle
    d      : screen distance
    """

    n = npx+1
    I = np.arange(1, n)
    x = I - n / 2
    y = n / 2 - I
    lam = 2.0*np.pi/k  # wavelength
    X = x[:, np.newaxis]
    Y = y[np.newaxis, :]

    ## diffraction mask
    M = np.zeros((n-1, n-1))  

    for i in range(nholes):
        t = i*2.0*np.pi/nholes
        xi = R*np.cos(t) 
        yi = R*np.sin(t) 
        ri = np.sqrt((X-xi)**2 + (Y-yi)**2 + d**2)
        M = M + (-1j/lam*np.exp(1j*k*ri)/ri)

    # (x, y) shifts by (u, v)
    xy = np.exp(1j*(u*X+v*Y))

    D1 = fftpack.fft2(xy*M)
    D2 = fftpack.fftshift(D1)

    abs_image = np.abs(D2)
    
    return abs_image

class Example(QtGui.QWidget):
    
    def __init__(self):
        super(Example, self).__init__()
        
        self.initUI()
        
    def initUI(self):
        
        kl = QtGui.QLabel(u'wavenumber: ')
        ul = QtGui.QLabel(u'x shift: ')
        vl = QtGui.QLabel(u'y shift: ')
        hl = QtGui.QLabel('pinholes: ')
        cl = QtGui.QLabel('palette: ')
        dl = QtGui.QLabel(u'screen distance:')
        rl = QtGui.QLabel(u'circle radius:')

        self.sk = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.sk.setToolTip(u'light <b>wavenumber</b>')
        self.sk.setRange(1,1000)
        self.sk.setValue(350)
        self.sk.valueChanged[int].connect(self.changeValue_k)
        self.kvl = QtGui.QLabel('350')

        self.su = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.su.setToolTip(u'<b>x</b> shift(%)')
        self.su.setRange(0,100)
        self.su.setValue(16)
        self.su.valueChanged[int].connect(self.changeValue_u)
        self.uvl = QtGui.QLabel('16 %')

        self.sv = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.sv.setToolTip(u'<b>y</b> shift(%)')
        self.sv.setRange(0,100)
        self.sv.setValue(16)
        self.sv.valueChanged[int].connect(self.changeValue_v)
        self.vvl = QtGui.QLabel('16 %')

        self.sd = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.sd.setToolTip(u'screen <b>distance</b>')
        self.sd.setRange(1,1000)  ## has to be divided by 100 to get d
        self.sd.setValue(10)
        self.sd.valueChanged[int].connect(self.changeValue_d)
        self.dvl = QtGui.QLabel('0.1')

        self.sr = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.sr.setToolTip(u'circle <b>radius</b>')
        self.sr.setRange(1,100)  ## has to be divided by 10 to get R
        self.sr.setValue(20)
        self.sr.valueChanged[int].connect(self.changeValue_r)
        self.rvl = QtGui.QLabel('2.0')

        self.hc = QtGui.QSpinBox(self)
        self.hc.setRange(1,25)
        self.hc.setValue(8)
        self.hc.setToolTip(u'number of diffraction pinholes')
        self.hc.valueChanged[int].connect(self.changeValue_holes)

        ## create palette
        self.cm = pg.GradientWidget(self, orientation='top')
        self.cm.setToolTip(u'right mouse click on strip: preset palettes\nleft mouse click: define new color\n'\
                            u'right mouse click on triangle: delete color\nleft mouse click on triangle: set color')
        self.cm.loadPreset('thermal')
        self.cm.setMinimumWidth(250)
        self.lut = self.cm.getLookupTable(100)
        self.cm.sigGradientChanged.connect(self.changePalette)

        ## create preview image
        self.p100 = pindr(npx=100)
        self.img100 = pg.ImageItem(self.p100)
        self.img100.setLookupTable(self.lut)
        pl = QtGui.QLabel(u'preview: ')
        pl.setToolTip(u'final colors may differ from these')
        self.pvw = pg.GraphicsView()
        self.pvw.setToolTip(u'final colors may slightly differ')
        self.pvw.addItem(self.img100)
        self.pvw.setMinimumSize(100,100)
        self.pvw.setMaximumSize(100,100)

        plot = QtGui.QPushButton('&Plot it')
        plot.clicked.connect(self.plotIt)
        save = QtGui.QPushButton('&Save')
        save.clicked.connect(self.showDialog)

        ## Create window with GraphicsView widget
        self.view = pg.GraphicsView()

        ## Create image item from numpy array
        self.p600 = pindr()
        self.img600 = pg.ImageItem(self.p600)
        self.img600.setLookupTable(self.lut)
        self.view.addItem(self.img600)
        self.view.scaleToImage(self.img600)
        self.view.setMinimumSize(600,600)
        self.view.setMaximumSize(600,600)

        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(kl, 0, 0)
        grid.addWidget(self.sk, 0, 1)
        grid.addWidget(self.kvl, 1, 1)

        grid.addWidget(ul, 2, 0)
        grid.addWidget(self.su, 2, 1)
        grid.addWidget(self.uvl, 3, 1)

        grid.addWidget(vl, 4, 0)
        grid.addWidget(self.sv, 4, 1)
        grid.addWidget(self.vvl, 5, 1)

        grid.addWidget(dl, 6, 0)
        grid.addWidget(self.sd, 6, 1)
        grid.addWidget(self.dvl, 7, 1)

        grid.addWidget(rl, 8, 0)
        grid.addWidget(self.sr, 8, 1)
        grid.addWidget(self.rvl, 9, 1)

        grid.addWidget(hl, 10, 0)
        grid.addWidget(self.hc, 10, 1)

        grid.addWidget(pl, 11, 0)
        ## grid.setRowStretch(11, 10)  # doesn't work
        grid.addWidget(self.pvw, 11, 1, alignment=QtCore.Qt.AlignHCenter)

        grid.addWidget(cl, 12, 0)
        grid.addWidget(self.cm, 12, 1)
        grid.addWidget(plot, 13, 1)
        grid.addWidget(save, 13, 0)


        grid.addWidget(self.view, 0, 2, 14, 1)
       
        self.setLayout(grid) 
        
        self.setGeometry(300, 300, 450, 400)
        self.setWindowTitle(u'Diffraction patterns')    
        self.show()

    def changeValue_k(self, value):
        self.p100 = pindr(k=value, u=self.su.value()*np.pi/50.0, v=self.sv.value()*np.pi/50.0, npx=100,
                     nholes=self.hc.value(), d=self.sd.value()/100.0, R=self.sr.value()/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)
        self.kvl.setText(str(value))

    def changeValue_u(self, value):
        self.p100 = pindr(k=self.sk.value(), u=value*np.pi/50.0, v=self.sv.value()*np.pi/50.0, npx=100,
                     nholes=self.hc.value(), d=self.sd.value()/100.0, R=self.sr.value()/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)
        self.uvl.setText(str(value)+" %")

    def changeValue_v(self, value):
        self.p100 = pindr(k=self.sk.value(), u=self.su.value()*np.pi/50.0, v=value*np.pi/50.0, npx=100,
                     nholes=self.hc.value(), d=self.sd.value()/100.0, R=self.sr.value()/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)
        self.vvl.setText(str(value)+" %")

    def changeValue_d(self, value):
        self.p100 = pindr(k=self.sk.value(), u=self.su.value()*np.pi/50.0, v=self.su.value()*np.pi/50.0, npx=100,
                     nholes=self.hc.value(), d=value/100.0, R=self.sr.value()/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)
        self.dvl.setText(str(value/100.0))

    def changeValue_r(self, value):
        self.p100 = pindr(k=self.sk.value(), u=self.su.value()*np.pi/50.0, v=self.su.value()*np.pi/50.0, npx=100,
                     nholes=self.hc.value(), d=value/100.0, R=value/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)
        self.rvl.setText(str(value/10.0))

    def changeValue_holes(self, value):
        self.p100 = pindr(k=self.sk.value(), u=self.su.value()*np.pi/50.0, v=self.sv.value()*np.pi/50.0, npx=100,
                     nholes=value, d=self.sd.value()/100.0, R=self.sr.value()/10.0)
        self.pvw.removeItem(self.img100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.addItem(self.img100)

    def plotIt(self):
        self.p600 = pindr(k=self.sk.value(), u=self.su.value()*np.pi/50.0, v=self.sv.value()*np.pi/50.0, npx=600,
                     nholes=self.hc.value(), d=self.sd.value()/100.0, R=self.sr.value()/10.0)
        ## rescale to the same LUT as preview image
        #a1 = np.amin(self.p100);  b1 = np.amax(self.p100)
        #a2 = np.amin(self.p600);  b2 = np.amax(self.p600)
        #c = (a1-b1)/(a2-b2)
        #d = a1 - c*a2
        #self.p600 = c*self.p600 + d
        ## even if you do the above, the two images will remain different for certain cases 
        ## -> linear rescaling doesn't help - histogram is different
        #print a1, b1
        #print np.amin(self.p600), np.amax(self.p600)
        self.view.removeItem(self.img600)
        self.img600.setImage(self.p600, lut=self.lut)
        self.view.scaleToImage(self.img600)
        self.view.addItem(self.img600)

    def changePalette(self):
        self.pvw.removeItem(self.img100)
        self.lut = self.cm.getLookupTable(100)
        self.img100.setImage(self.p100, lut=self.lut)
        self.pvw.scaleToImage(self.img100)
        self.pvw.addItem(self.img100)

    def showDialog(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, u'Save image', 'image.png')
        self.img600.save(fname)


        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
