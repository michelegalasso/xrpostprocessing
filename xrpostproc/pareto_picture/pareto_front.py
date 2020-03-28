"""
@file:      pareto_front.py
@author:    Michele Galasso
@contact:   michele.galasso@skoltech.ru
@date:      25 March 2020
@brief:     Function for finding a 2D Pareto frontier given two lists of matched length.
"""


def pareto_front(Xs, Ys, maxX=True, maxY=True):
    """
    Function for finding a 2D Pareto frontier given two lists of matched length.
    :Link: http://code.activestate.com/recipes/578230-pareto-front/
    """
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]
    rest = []
    for pair in myList[1:]:
        if maxY:
            if pair[1] >= p_front[-1][1]:
                p_front.append(pair)
            else:
                rest.append(pair)
        else:
            if pair[1] <= p_front[-1][1]:
                p_front.append(pair)
            else:
                rest.append(pair)
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    restX = [pair[0] for pair in rest]
    restY = [pair[1] for pair in rest]
    return p_frontX, p_frontY, restX, restY
