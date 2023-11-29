#!/home/xihe/.conda/envs/physics/bin/python

import uproot
import plotext as plt
import argparse

def plot_histogram(filename, hist_name=None):
    # Load the .root file
    file = uproot.open(filename)
    
    # If no histogram name is provided, default to the first histogram
    if hist_name is None:
        # Get the first histogram key
        for key in file.keys():
            if key.endswith(';1'):
                hist_name = key
                break
        else:
            raise ValueError("No histograms found in the file")
    
    # Access the histogram
    histogram = file[hist_name]

    
    # Get the histogram data
    values, edges = histogram.to_numpy()
    
    # Get the bin centers
    bin_centers = (edges[:-1] + edges[1:]) / 2
    
    # Plot the histogram using plotext
    plt.plot(bin_centers, values)
    
    # Set the title and axis labels from the ROOT histogram
    plt.title(histogram.title)
    plt.xlabel((histogram.axis(0).labels()))
    plt.ylabel((histogram.axis(-1).labels()))
    
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot a histogram from a ROOT file.')
    parser.add_argument('filename', help='The name of the ROOT file.')
    parser.add_argument('--hist_name', default=None, help='The name of the histogram to plot. If omitted, the first histogram will be plotted.')
    
    args = parser.parse_args()
    
    plot_histogram(args.filename, args.hist_name)
