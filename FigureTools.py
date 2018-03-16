from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib.pyplot as plt

def figure_numbering(text): 
    artist = AnchoredText(text,
                      prop=dict(size=12), 
                      frameon=True,
                      loc=2, #  2: upper left 
                      pad=0.1,
                      bbox_to_anchor=(-.1, 1.15),
                      bbox_transform=plt.gca().transAxes
                      )
    #artist.patch.set_boxstyle(fc='#FFFFFF00',ec='#FFFFFF00')
    artist.patch.set_fc('white')
    artist.patch.set_ec('white')
    return artist