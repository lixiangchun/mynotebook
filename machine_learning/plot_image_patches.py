import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot(x, ksize_rows=299, ksize_cols=299):
    nr = x.shape[1]
    nc = x.shape[2]
    fig = plt.figure(figsize=(nr, nc))
    gs = gridspec.GridSpec(nr, nc)
    gs.update(wspace=0.05, hspace=0.05)

    for i in range(nr):
        for j in range(nc):
            ax = plt.subplot(gs[i*nr+j])
            plt.axis('off')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            #ax.set_aspect('equal')
            #plt.imshow(sample.reshape(ksize_rows, ksize_cols), cmap='Greys_r')
            plt.imshow(x[0,i,j,].reshape(ksize_rows, ksize_cols, 3))

    return fig

#fig = plot(x)
#plt.savefig('out/{}.png'.format(str(i).zfill(3)), bbox_inches='tight')
#plt.close(fig)
