import numpy as np
from functools import partial
import math

def trapezoid_equidistant(function,bounds,h=30):
    """
    Izračuna določeni integral podane funkcije med dvema mejama z uporabo trapeznega pravila.

    Parameters:
        function (callable): funkcija, ki se integrira.
        bounds (tuple): integracijske meje (nižja, višja).
        h (int): število ekvidistantnih točk, privzeto 30.

    Returns:
        float: izračunani približek vrednosti integrala.
    """
    #Calculate points
    points, distance = np.linspace(bounds[0],bounds[1],h,retstep=True)
    #Return sum
    return distance/2 * np.sum([function(points[i])+function(points[i+1]) for i in range(len(points)-1)])

def bessel_0(x):
    """
   Izračuna vrednost besselove funkcije ničtega reda.

   Parameters:
       x (float): vhodna vrednost.

   Returns:
       float: vrednost besselove funkcije ničtega reda.
   """
    return trapezoid_equidistant(partial(internal_function,v=x),(0,math.pi),h=1000)

def internal_function(x,v):
    """
    Funkcija, ekvivalentna notranji funkciji besselovega integrala.

    Parameters:
        x (float): parameter t - iz besselove funkcije.
        v (float): parameter x - iz besselove funkcije.

    Returns:
        float: Vrednost funkcije.
    """
    return np.cos(v * np.sin(x))/(math.pi)


def point_from_bezier(control_points, t):
    """
    Izračuna vrednost točke na Bezierjevi krivulji, podane s kontrolnimi točkami, na mestu t.

    Funkcija rekurzivno izračuna točko na Bezierjevi krivulji z uporabo De Casteljau-jevega algoritma.

    Parameters:
       control_points (list): list parov števil, ki predstavljajo kontrolne točke (x, y).
       t (float): parametrična vrednost med 0 in 1.

    Returns:
       tuple: par števil (x, y), ki predstavlja točko na krivulji za parameter p.
    """
    if (len(control_points) == 1):
        return control_points[0]

    cpoints_new = [
        (
            (1 - t) * control_points[i][0] + t * control_points[i + 1][0],
            (1 - t) * control_points[i][1] + t * control_points[i + 1][1]
        )
        for i in range(len(control_points) - 1)
    ]

    return point_from_bezier(cpoints_new, t)

def shoelace_area(points):
    """
    Izračuna ploščino poligona, omejenega z točkami, z uporabo algoritma 'Shoelace'.

    Parameters:
        points (list): točke, ki omejujejo poligon.

    Returns:
        float: izračunano ploščina poligona.
    """
    n = len(points)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += (points[i][0] * points[j][1]) - (points[j][0] * points[i][1])
    return area/2