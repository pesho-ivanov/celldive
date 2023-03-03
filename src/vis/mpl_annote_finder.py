# from http://scipy-cookbook.readthedocs.io/items/Matplotlib_Interactive_Plotting.html

import math
import matplotlib.pyplot as plt

class AnnoteFinder(object):
    """callback for matplotlib to display an annotation when points are
    clicked on.  The point which is closest to the click and within
    xtol and ytol is identified.
    
    Register this function like this:
    
    scatter(xdata, ydata)
    af = AnnoteFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self, xdata, ydata, annotes, ax=None, xtol=None, ytol=None, initOn=False, someOn=None):
        self.data = list(zip(xdata, ydata, annotes))
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if ax is None: self.ax = plt.gca()
        else: self.ax = ax
        self.drawnAnnotations = {}
        self.links = []
        if initOn:
            for x, y, annote in self.data:
                self.drawAnnote(ax, x, y, annote)
        if someOn:
            #data = [e + (z,) for (e, z) in zip(x, y)]
            ext_data = list(zip(xdata, ydata, annotes, someOn))
            for x, y, annote, on in ext_data:
                if on:
                    self.drawAnnote(ax, x, y, annote)

    def distance(self, x1, x2, y1, y2):
        return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2))

    def __call__(self, event):
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            if (self.ax is None) or (self.ax is event.inaxes):
                annotes = []
                # print(event.xdata, event.ydata)
                for x, y, a in self.data:
                    # print(x, y, a)
                    if ((clickX-self.xtol < x < clickX+self.xtol) and
                            (clickY-self.ytol < y < clickY+self.ytol)):
                        annotes.append(
                            (self.distance(x, clickX, y, clickY), x, y, a))
                if annotes:
                    annotes.sort()
                    distance, x, y, annote = annotes[0]
                    self.drawAnnote(event.inaxes, x, y, annote)
                    for l in self.links:
                        l.drawSpecificAnnote(annote)

    def drawAnnote(self, ax, x, y, annote):
        """
        Draw the annotation on the plot
        """
        if (x, y) in self.drawnAnnotations:
            markers = self.drawnAnnotations[(x, y)]
            for m in markers:
                #print('bla')
                m.set_visible(not m.get_visible())
            self.ax.figure.canvas.draw_idle()
        else:
            if 'TRB' in annote:
                ha = 'left'
            else:
                ha = 'right'
            t = ax.text(x, y, " %s" % (annote), size=7, ha=ha)
            #m = ax.scatter([x], [y], marker='d', c='r', zorder=100)
            #self.drawnAnnotations[(x, y)] = (t, m)
            self.drawnAnnotations[(x, y)] = (t,)
            self.ax.figure.canvas.draw_idle()

    def drawSpecificAnnote(self, annote):
        for x, y, a in self.data:
            if a == annote:
                self.drawAnnote(self.ax, x, y, a)
            
    @staticmethod
    def linkAnnotationFinders(afs):
        '''Usage: linkAnnotationFinders([af1, af2])'''
        for i in range(len(afs)):
            allButSelfAfs = afs[:i]+afs[i+1:]
            afs[i].links.extend(allButSelfAfs)
            
# usage
def demo():
    x = range(10)
    y = range(10)
    annotes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']

    fig = plt.figure()
    
    # first subplot
    ax = fig.add_subplot(121)
    ax.scatter(x,y)
    af1 = AnnoteFinder(x, y, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af1)
    
    # second subplot
    ax = fig.add_subplot(122)
    ax.scatter(x,y)
    af2 = AnnoteFinder(x, y, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af2)
    
    AnnoteFinder.linkAnnotationFinders([af1, af2])
    
    plt.show()