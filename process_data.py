# Importing the libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.widgets import Slider
from scipy.optimize import leastsq
import os 

# Plots settings 
plt.rc("font", size=10)
plt.rcParams["font.family"] = "Times New Roman"

def coords_reproject(cluster, coords, gaia_mag, plots=True, members=False):
    # Importing gaia and uv data
    print(f"\nImporting data for {cluster}...\n")
    gaia = pd.read_csv(f'input/star_table_{cluster}.csv')
    uv = pd.read_csv(f'input/uv.csv')

    if members:
        # Importing members data from Hunt catalog 
        print(f"\nImporting members data...")
        members_m = pd.read_csv('../MS_members/members_m_final.csv')
        members_p = pd.read_csv('../MS_members/members_p_final.csv')
        print(f"Members data imported.\n")

        members_comb = pd.concat([members_m, members_p], ignore_index=True)

        members_comb_filtered = members_comb[members_comb['name'] == cluster]

        if members_comb_filtered.empty:
            print("Error: No members data found for the given cluster. Try to look for cluster name in Hunt catalogue. Continue without plotting members.")
            members = False

    # Adding id column to gaia and uv data (from 1 to len(data))
    gaia['id'] = range(1, len(gaia) + 1)   
    uv['id'] = range(1, len(uv) + 1)   

    print(f"\nData for {cluster} imported.")
    print(f"Number of stars in gaia data: {len(gaia)}")
    print(f"Number of stars in uv data: {len(uv)}\n")

    # Selecting columns for gaia and uv data and saving them to output folder
    print(f"\nSaving selected columns from gaia and uv data to output folder...")
    gaia_selected = gaia[['id', 'ra', 'dec', gaia_mag]]
    gaia_selected.to_csv(f'output/coords/gaia_sky.dat', sep='\t', index=False, header=False)

    uv_selected = uv[['id', 'ra', 'dec', 'u']]
    uv_selected.to_csv(f'output/coords/uv_sky.dat', sep='\t', index=False, header=False)
    print(f"Selected columns saved to output folder.\n")

    # Changing coordinates from RA-DEC to x-y
    print(f"\nChanging coordinates from RA-DEC to x-y...")
    os.system(f'project_coords output/coords/gaia_sky.dat 1 2 {coords[0]} {coords[1]} asec outfile=output/coords/gaia_cart.dat')
    os.system(f'project_coords output/coords/uv_sky.dat 1 2 {coords[0]} {coords[1]} asec outfile=output/coords/uv_cart.dat')
    print(f"Coordinates changed from RA-DEC to x-y.\n")

    # Importing datasets with cartesian coordinates
    print(f"\nImporting datasets with cartesian coordinates...")
    gaia_cart = pd.read_csv(f'output/coords/gaia_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', gaia_mag])
    uv_cart = pd.read_csv(f'output/coords/uv_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])
    print(f"Datasets with new coordinates imported.\n")

    # Plotting RA-DEC and X-Y coords
    if plots:
        print(f"\nPlotting RA-DEC and X-Y coords...")
        plot_sky(gaia, uv, members_comb_filtered, members)
        plot_cart(gaia_cart, uv_cart)
        print(f"Plots saved to plots folder.\n")

def matching(cluster, gaia_mag, matchrad=3.0, trirad=0.001, nobj=40, plots=True):
    # Matching the datasets
    #http://spiff.rit.edu/match/match-0.16/match.html

    # parameters like: matchrad=3.0 trirad=0.001 nobj=40 recalc is changable and requares playing with them to get the best result
    print(f"\nMatching stars from Gaia and UV datasets for {cluster}...\n")
    os.system(f'match output/coords/gaia_cart.dat 1 2 3 output/coords/uv_cart.dat 1 2 3 outfile=output/matched/matched id1=0 id2=0 matchrad={matchrad} trirad={trirad} nobj={nobj} recalc')
    print(f"\nStars matched.\n")

    # Importing matched datasets
    print(f"\nImporting matched datasets for {cluster}...")
    gaia_matched = pd.read_csv(f'output/matched/matched.mtA', sep='\s+', header=None, names=['id', 'x', 'y', gaia_mag]) 

    uv_matched = pd.read_csv(f'output/matched/matched.mtB', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])
    uv_unmatched = pd.read_csv(f'output/matched/matched.unB', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])
    print(f"Matched datasets imported.\n")

    # Importing original gaia data for sining SOURCE_IDs to matched data
    print(f"\nImporting original gaia data for {cluster}...")
    gaia = pd.read_csv(f'input/star_table_{cluster}.csv')
    # Adding id column to gaia and uv data (from 1 to len(data))
    gaia['id'] = range(1, len(gaia) + 1)  
    print(f"Original gaia data imported.\n")

    # Importing matched datasets 
    columns_to_add = ['SOURCE_ID', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'ra', 'dec', 'parallax'] # List of columns to add
    gaia_matched = gaia_matched.drop(columns=gaia_mag) # Remove the column represented by the variable gaia_mag
    gaia_matched = pd.merge(gaia_matched, gaia[['id'] + columns_to_add], on='id') # Merge the dataframes on the 'id' column
    gaia_matched['SOURCE_ID'] = gaia_matched['SOURCE_ID'].astype('int64')

    uv_matched = pd.concat([uv_matched, gaia_matched[columns_to_add]], axis=1)
    uv_matched['SOURCE_ID'] = uv_matched['SOURCE_ID'].astype('int64')

    # Saving matched datasets
    print(f"\nSaving matched datasets for {cluster}...")
    gaia_matched.to_csv(f'output/matched/gaia_matched.dat', sep='\t', index=False)
    uv_matched.to_csv(f'output/matched/uv_matched.dat', sep='\t', index=False)
    print(f"Matched datasets saved.\n")

    # Plotting X-Y positions of matched stars
    if plots:
        # Importing datasets with cartesian coordinates
        print(f"\nImporting datasets with cartesian coordinates...")
        gaia_cart = pd.read_csv(f'output/coords/gaia_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', gaia_mag])
        uv_cart = pd.read_csv(f'output/coords/uv_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])
        print(f"Datasets with new coordinates imported.\n")

        print(f"\nPlotting matched stars for {cluster}...")
        plot_matched(gaia_matched, uv_matched, uv_unmatched, gaia_cart, uv_cart)

def filtering(cluster, plots=True):
    # Importing members data from Hunt catalog
    print(f"\nImporting members data...")
    members_m = pd.read_csv('../MS_members/members_m_final.csv')
    members_p = pd.read_csv('../MS_members/members_p_final.csv')
    print(f"Members data imported.\n")

    # Importing matched datasets
    print(f"\nImporting matched datasets for {cluster}...")
    gaia_matched = pd.read_csv(f'output/matched/gaia_matched.dat', sep='\t')
    uv_matched = pd.read_csv(f'output/matched/uv_matched.dat', sep='\t')

    # Filtering members
    members_comb = pd.concat([members_m, members_p], ignore_index=True)
    members_comb_filtered = members_comb[members_comb['probability'] >= 0.7]

    # Filtering gaia data
    mask = gaia_matched['SOURCE_ID'].isin(members_comb_filtered['SOURCE_ID'])
    gaia_filtered = gaia_matched[mask]
    gaia_filtered['SOURCE_ID'] = gaia_filtered['SOURCE_ID'].astype('int64')
    
    # Filtering uv data
    mask = uv_matched['SOURCE_ID'].isin(members_comb_filtered['SOURCE_ID'])
    uv_filtered = uv_matched[mask]
    uv_filtered['SOURCE_ID'] = uv_filtered['SOURCE_ID'].astype('int64')

    # Saving filtered datasets
    print(f"\nSaving filtered datasets for {cluster}...")
    gaia_filtered.to_csv(f'output/matched/gaia_filtered.dat', sep='\t', index=False)
    uv_filtered.to_csv(f'output/matched/uv_filtered.dat', sep='\t', index=False)
    print(f"Filtered datasets saved.\n")

    # Plotting X-Y positions of filtered stars
    if plots: 
        print(f"\nImporting datasets with cartesian coordinates...")
        gaia_cart = pd.read_csv(f'output/coords/gaia_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', 'phot_g_mean_mag'])
        uv_cart = pd.read_csv(f'output/coords/uv_cart.dat', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])
        print(f"Datasets with new coordinates imported.\n")
         
        print(f"\nPlotting filtered stars for {cluster}...")
        plot_filtered(gaia_filtered, uv_filtered, gaia_cart, uv_cart)

def diagram(cluster, colour1, colour2, adjust=False, cmd=False):
    # Importing main sequence data
    print(f"\nImporting main sequence data...")
    ms1 = pd.read_csv('../MS_members/lines_UBPRPG_1.txt', delimiter = '\s+')
    ms2 = pd.read_csv('../MS_members/lines_UBPRPG_2.txt', delimiter = '\s+')
    print(f"Main sequence data imported.\n")

    if cmd:
        pass
    else:
        if colour1 in ms1.columns and colour2 in ms1.columns:
            ms = ms1[[colour1, colour2]]
        elif colour1 in ms2.columns and colour2 in ms2.columns:
            ms = ms2[[colour1, colour2]]
        else:
            print(f"Combination of colours {colour1} and {colour2} cant be used, because filters are located in different files. Try another combination.")

        # Fitting MS
        coeffs = np.polyfit(ms[colour1], ms[colour2], 9)
        poly = np.poly1d(coeffs)
        ms_fit = poly(ms[colour1])

    # Importing filtered datasets
    print(f"\nImporting filtered datasets for {cluster}...")
    gaia_filtered = pd.read_csv(f'output/matched/gaia_filtered.dat', sep='\t')
    uv_filtered = pd.read_csv(f'output/matched/uv_filtered.dat', sep='\t')
    print(f"Filtered datasets imported.\n")

    # Importing matched datasets
    print(f"\nImporting matched datasets for {cluster}...")
    gaia_matched = pd.read_csv(f'output/matched/gaia_matched.dat', sep='\t')
    uv_matched = pd.read_csv(f'output/matched/uv_matched.dat', sep='\t')
    print(f"Matched datasets imported.\n")

    # Calculating colors
    print(f"\nCalculating for {colour1} and {colour2}...")
    gaia_filter_1 = f'phot_{colour1.split("-")[0].lower()}_mean_mag'
    gaia_filter_2 = f'phot_{colour1.split("-")[1].lower()}_mean_mag'

    gaia_filtered[colour1] = gaia_filtered[gaia_filter_1] - gaia_filtered[gaia_filter_2]
    gaia_matched[colour1] = gaia_matched[gaia_filter_1] - gaia_matched[gaia_filter_2]

    if cmd: 
        if colour2 == 'U':
            uv_filtered[colour2] = uv_filtered['u'] 
            uv_matched[colour2] = uv_matched['u'] 
        else:
            gaia_filter_3 = f'phot_{colour2.lower()}_mean_mag'

            gaia_filtered[colour2] = gaia_filtered[gaia_filter_3]
            gaia_matched[colour2] = gaia_matched[gaia_filter_3]
    else:
        gaia_filter_3 = f'phot_{colour2.split("-")[1].lower()}_mean_mag'

        uv_filtered[colour2] = uv_filtered['u'] - gaia_filtered[gaia_filter_3]
        uv_matched[colour2] = uv_matched['u'] - gaia_matched[gaia_filter_3]

    print(f"Colours calculated.\n")

    # Plotting colour-colour diagram
    if cmd: 
        print(f"\nPlotting colour-magnitude diagram for {cluster} with colours {colour1} and {colour2}...")
        if not os.path.isdir(f'plots/CMD'): os.system(f'mkdir -p plots/CMD')
        plot_cmd(gaia_matched, uv_matched, gaia_filtered, uv_filtered, colour1, colour2)
        print(f"Colour-magnitude diagram saved.\n")
    else:
        print(f"\nPlotting colour-colour diagram for {cluster} with colours {colour1} and {colour2}...")
        if not os.path.isdir(f'plots/CC'): os.system(f'mkdir -p plots/CC')
        plot_cc(gaia_matched, uv_matched, gaia_filtered, uv_filtered, ms, colour1, colour2, adjust)
        print(f"Colour-colour diagram saved.\n")

def radius(ra, dec):
    ra_min = np.min(ra)
    ra_max = np.max(ra)
    radius_ra = (ra_max - ra_min)/2

    dec_min = np.min(dec)
    dec_max = np.max(dec)
    radius_dec = (dec_max - dec_min)/2

    return [radius_ra, radius_dec]

def center(ra, dec):
    mean_ra = ra.mean()
    mean_dec = dec.mean()
    return [mean_ra, mean_dec]

def plot_sky(gaia, uv, members_comb_filtered, members=False):
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel(r"$RA$ [deg]")
    ax.xaxis.label.set_fontsize(25)
    ax.set_ylabel(r"$DEC$ [deg]")
    ax.yaxis.label.set_fontsize(25)
    ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
    plt.tight_layout()  

    ax.scatter(gaia['ra'], gaia['dec'], color='red', s=10, alpha=0.5, label='Gaia (sky)')
    ax.scatter(uv['ra'], uv['dec'], color='blue', s=10, alpha=0.5, label='UV (sky)')

    if members:
        ax.scatter(members_comb_filtered['ra'], members_comb_filtered['dec'], color='black', s=20, marker='^', label='Members')

    ax.add_patch(Ellipse((center(gaia['ra'], gaia['dec'])), 2*radius(gaia['ra'], gaia['dec'])[0], 2*radius(gaia['ra'], gaia['dec'])[1], edgecolor='red', fc='None', lw=2))
    ax.add_patch(Ellipse((center(uv['ra'], uv['dec'])), 2*radius(uv['ra'], uv['dec'])[0], 2*radius(uv['ra'], uv['dec'])[1], edgecolor='blue', fc='None', lw=2))

    ax.legend(loc='upper right', fontsize=15)

    fig.savefig(f'plots/sky_coords.png', bbox_inches='tight')

def plot_cart(gaia_cart, uv_cart):
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel(r"$X$ [arcsec]")
    ax.xaxis.label.set_fontsize(25)
    ax.set_ylabel(r"$Y$ [arcsec]")
    ax.yaxis.label.set_fontsize(25)
    ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
    plt.gca().set_aspect('equal')
    plt.tight_layout()  

    ax.scatter(gaia_cart['x'], gaia_cart['y'], color='red', s=10, alpha=0.5, label='Gaia (cart)')
    ax.scatter(uv_cart['x'], uv_cart['y'], color='blue', s=10, label='UV (cart)')

    ax.legend(loc='upper right', fontsize=15)

    fig.savefig(f'plots/cart_coords.png', bbox_inches='tight')

def plot_matched(gaia_matched, uv_matched, uv_unmatched, gaia_cart, uv_cart):
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel(r"$X$ [arcsec]")
    ax.xaxis.label.set_fontsize(25)
    ax.set_ylabel(r"$Y$ [arcsec]")
    ax.yaxis.label.set_fontsize(25)
    ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
    plt.gca().set_aspect('equal')
    plt.tight_layout()  

    ax.scatter(gaia_cart['x'], gaia_cart['y'], color='red', s=5, label=f'Gaia {len(gaia_cart)}', alpha=0.3)
    ax.scatter(uv_cart['x'], uv_cart['y'], color='blue', s=5, label=f'UV {len(uv_cart)}', alpha=0.3)

    ax.scatter(gaia_matched['x'], gaia_matched['y'], color='red', s=10, marker='^', label=f'Gaia (matched) {len(gaia_matched)}')
    ax.scatter(uv_matched['x'], uv_matched['y'], color='blue', s=10, marker='^', label=f'UV (matched) {len(uv_matched)}')

    ax.scatter(uv_unmatched['x'], uv_unmatched['y'], color='green', s=10, marker='x', label=f'UV (unmatched) {len(uv_unmatched)}')

    ax.legend(loc='upper right', fontsize=15)
    
    fig.savefig(f'plots/cart_coords_matched.png', bbox_inches='tight')

def plot_filtered(gaia_filtered, uv_filtered, gaia_cart, uv_cart):
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel(r"$X$ [arcsec]")
    ax.xaxis.label.set_fontsize(25)
    ax.set_ylabel(r"$Y$ [arcsec]")
    ax.yaxis.label.set_fontsize(25)
    ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
    plt.tight_layout()  

    ax.scatter(gaia_cart['x'], gaia_cart['y'], color='red', s=5, label=f'Gaia {len(gaia_cart)}', alpha=0.2)
    ax.scatter(uv_cart['x'], uv_cart['y'], color='blue', s=5, label=f'UV {len(uv_cart)}', alpha=0.3)

    ax.scatter(gaia_filtered['x'], gaia_filtered['y'], color='red', s=20, marker='^', label=f'Gaia (filtered) {len(gaia_filtered)}')
    ax.scatter(uv_filtered['x'], uv_filtered['y'], color='blue', s=20, marker='^', label=f'UV (filtered) {len(uv_filtered)}')

    ax.legend(loc='upper right', fontsize=15)

    fig.savefig(f'plots/cart_coords_filtered.png', bbox_inches='tight')

def plot_cc(gaia_matched, uv_matched, gaia_filtered, uv_filtered, ms, colour1, colour2, adjust):
    fig, ax = plt.subplots(figsize=(8, 8))

    if adjust: 
        def update(val):
            x = sx.val
            y = sy.val
            scatter_ms.set_offsets(np.c_[ms[colour1] + x, ms[colour2] + y])
            ax.legend(loc='upper right', fontsize=15)
            fig.canvas.draw_idle()
            
        plt.subplots_adjust(left=0.1, bottom=0.25)

        # Set labels and other aesthetics
        ax.set_xlabel(f"{colour1} [mag]")
        ax.set_ylabel(f"{colour2} [mag]")
        ax.invert_yaxis()
        ax.tick_params(axis="both", which="major", labelsize=20)
        ax.tick_params(axis="both", which="minor", labelsize=20)
        plt.grid(alpha=0.3)

        # Initial data points
        scatter_matched = ax.scatter(gaia_matched[colour1], uv_matched[colour2], color='green', s=7, alpha=0.3, label='Matched stars')
        scatter_members = ax.scatter(gaia_filtered[colour1], uv_filtered[colour2], color='magenta', s=17, marker='^', label=f'Membership')
        scatter_ms = ax.scatter(ms[colour1], ms[colour2], color='red', s=5, label='Main sequence')

        # Add sliders for interactive update
        axcolor = 'lightgoldenrodyellow'
        ax_x = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
        ax_y = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

        sx = Slider(ax_x, 'X', -5, 5, valinit=0)
        sy = Slider(ax_y, 'Y', -5, 5, valinit=0)

        sx.on_changed(update)
        sy.on_changed(update)

        # Show plot
        plt.show()

        # Save plot
        x_new = sx.val
        y_new = sy.val

        fig2, ax2 = plt.subplots(figsize=(8, 8))
        ax2.set_xlabel(f"{colour1} [mag]")
        ax2.xaxis.label.set_fontsize(25)
        ax2.set_ylabel(f"{colour2} [mag]")
        ax2.yaxis.label.set_fontsize(25)
        ax2.invert_yaxis()
        ax2.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
        ax2.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
        plt.grid(alpha=0.3)
        plt.tight_layout()  

        ax2.scatter(gaia_matched[colour1], uv_matched[colour2], color='green', s=7, alpha=0.3, label='Matched stars')
        ax2.scatter(gaia_filtered [colour1], uv_filtered [colour2], color='magenta', s=17, marker='^', label='Membership')

        ax2.plot(ms[colour1]+x_new, ms[colour2]+y_new, color='red', linewidth=2, label=f'Main sequence ({colour1}={x_new:.2}, {colour2}={y_new:.2})')
        ax2.legend(loc='upper right', fontsize=15)

        fig2.savefig(f'plots/CC/{colour1}_{colour2}.png', bbox_inches='tight')
    else:
        ax.set_xlabel(f"{colour1} [mag]")
        ax.xaxis.label.set_fontsize(25)
        ax.set_ylabel(f"{colour2} [mag]")
        ax.yaxis.label.set_fontsize(25)
        ax.invert_yaxis()
        ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
        ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
        plt.grid(alpha=0.3)
        plt.tight_layout()  

        ax.scatter(gaia_matched[colour1], uv_matched[colour2], color='green', s=7, alpha=0.3, label='Matched stars')
        ax.scatter(gaia_filtered[colour1], uv_filtered[colour2], color='magenta', s=17, marker='^', label=f'Membership')

        ax.plot(ms[colour1], ms[colour2], color='red', linewidth=2, label='Main sequence')
        ax.legend(loc='upper right', fontsize=15)

        fig.savefig(f'plots/CC/{colour1}_{colour2}.png', bbox_inches='tight')

def plot_cmd(gaia_matched, uv_matched, gaia_filtered, uv_filtered, colour1, colour2):
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel(f"{colour1} [mag]")
    ax.xaxis.label.set_fontsize(25)
    ax.set_ylabel(f"{colour2} [mag]")
    ax.yaxis.label.set_fontsize(25)
    ax.invert_yaxis()
    ax.tick_params(axis="both", which="major", length=10, width=0.5, labelsize=20)
    ax.tick_params(axis="both", which="minor", length=5, width=0.5, labelsize=20)
    plt.grid(alpha=0.3)
    plt.tight_layout()  

    if colour2 == 'U':
        ax.scatter(gaia_matched[colour1], uv_matched[colour2], color='green', s=7, alpha=0.3, label='Matched stars')
        ax.scatter(gaia_filtered[colour1], uv_filtered[colour2], color='magenta', s=17, marker='^', label=f'Membership')
    else:
        ax.scatter(gaia_matched[colour1], gaia_matched[colour2], color='green', s=7, alpha=0.3, label='Matched stars')
        ax.scatter(gaia_filtered[colour1], gaia_filtered[colour2], color='magenta', s=17, marker='^', label=f'Membership')

    ax.legend(loc='upper right', fontsize=15)

    fig.savefig(f'plots/CMD/{colour1}_{colour2}.png', bbox_inches='tight')