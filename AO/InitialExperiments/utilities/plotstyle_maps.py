from matplotlib import colors

def get_color_map():
    """Specify colors to use for visualizing mnist"""
    cmap = {}
    cmap['color_black'] = '#5A5A5A'
    cmap['color_red'] = '#CD3333'
    cmap['color_yel'] = '#E3CF57'
    cmap['color_cyan'] = '#42c9c9ff'
    cmap['color_lightblue'] = '#3380f2'
    cmap['color_darkgray'] = '#838B8B'
    cmap['color_green'] = '#aaffdd' #'#aaffdd' # #003700
    cmap['color_lime'] = '#ddffbb'
    cmap['color_pink'] = '#FFAEB9'
    cmap['color_lightgreen'] = '#c0ff0c'
    cmap['color_orange'] = '#f5a565'
    cmap['color_darkgreen'] = '#40826d'
    cmap['color_darkblue'] = '#00008B'
    cmap['white'] = '#FFFFFF'
    cmap['transparent'] = '#FFFFFF00'
    return cmap

def get_dash_map():
    # dashes=(dash_len, spacing)
    dash_map = {}
    dash_map['dash1'] = (10, 0)
    dash_map['dash2'] = (2, 2)
    dash_map['dash3'] = (4, 2)
    dash_map['dash4']  = (3, 2, 1, 2, 1, 2)
    dash_map['dash5'] = (3, 1, 1, 1, 1, 1)
    dash_map['dash6'] = (10, 2, 2, 2, 2, 2)
    dash_map['dash7'] = (5, 1, 3, 1)
    dash_map['dash8'] = (10, 1)
    dash_map['dash9'] = (20, 2)
    return dash_map

def get_marker_map():
    marker_map = {}
    marker_map['marker1'] = 'o'
    marker_map['markerNone'] = None
    marker_map['marker2'] = 'o'
    marker_map['marker3'] = 's'
    marker_map['marker4']  = 'v'
    marker_map['marker5'] = 'p'
    marker_map['marker6'] = 'd'
    marker_map['marker7'] = 'h'
    marker_map['marker8'] = 'D'
    marker_map['marker9'] = 'P'
    return marker_map
