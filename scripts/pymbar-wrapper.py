#!/usr/bin/python
#
# Wrapper for PyMBAR
# expects metadata and data files in the Grossman WHAM format
#

import os
import sys
import optparse
from math import pi
import random
import numpy

from pymbar import timeseries # timeseries analysis
from pymbar import pymbar

# Constants.
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K
temperature = 300 # assume a single temperature -- can be overridden with data from center.dat 
# beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
N_max = 501 # maximum number of snapshots/simulation

def file_len(fname, ignore='#'):
    """ Count the number of lines in a file """
    with open(fname) as f:
        comments = 0
        for i, l in enumerate(f):
            if l.startswith(ignore):
                comments += 1
            pass
    return i + 1 - comments

def main():    
    usage = """
        usage: %prog [options] <metadata file>
    """
    
    parser = optparse.OptionParser(usage)
    #parser.add_option("-m", "--config", dest="config_file", default="config.ini", help="Config file [default: %default]")
    #parser.add_option("-p", "--percent", dest="percent", default=100, type="int", help="Percent of data to extract for WHAM [default: %default %%]")
    #parser.add_option("-r", "--random", dest="random", default=False, action="store_true", help="Select data points at random [default: %default]")
    #parser.add_option("-s", "--start-index", dest="start_index", default=0, type="int", help="Starting index for column data [default: %default]")
    
    parser.add_option("-o", "--outfile", dest="output_file", default='mbar_pmf.out', help="Output file for PMF [default: %default]")
    parser.add_option("-t", "--temperature", dest="temperature", default=300., type="float", help="Initial temperature in K [default: %default K]")
    
    #parser.add_option("--min", dest="min_coord", default=0, type="int", help="Min coordinate for 1D PMF [default: %default]")
    #parser.add_option("--max", dest="max_coord", default=360, type="int", help="Max coordinate for 1D PMF [default: %default]")
    parser.add_option("-b", "--bins", dest="bins", default=50, type="int", help="Number of bins for 1D PMF [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if args < 1:
        parser.error('No metadata file passed')
    elif not os.path.exists(args[0]):
        parser.error('Metadata file not found')
    
    metadata = [] # stores metadata per umbrella

    # open the wham metadata file
    print "Opening metadata file %s" % args[0]
    f = open(args[0], 'r')
    metadata_lines = f.readlines()
    f.close()
    
    metadata = []
    N_max = 0 # the max number of snapshots per umbrella
    different_temperatures = False
    
    # first get all the metadata and count the max number of snapshots per umbrella
    for line in metadata_lines:
        if line.startswith('#'):
            continue
        clean_split = filter(None, line.strip().expandtabs().split(' '))
        if not os.path.exists(clean_split[0]):
            print "Data file %s doesn't exist, skipping this replica" % clean_split[0]
            continue
        else:
            nsnapshots = file_len(clean_split[0])
            if nsnapshots > N_max:
                N_max = nsnapshots
            
            # /path/to/timeseries/file  loc_win_min spring  [correl time] [temperature]
            current_meta = { 'path': clean_split[0], 'coord': float(clean_split[1]), 'k': float(clean_split[2]), 'n': nsnapshots }
            # TODO: spring constant units?
            #   K_k[k] = float(tokens[1]) * (numpy.pi/180)**2 # spring constant (read in kJ/mol/rad**2, converted to kJ/mol/deg**2)    
            
            if len(clean_split) >= 4:
                # TODO: temperature the 4rd or 5th value???
                # temperature might be the 4th value...
                current_meta['t'] = float(clean_split[3])
                different_temperatures = True
            metadata.append(current_meta)
            
    print "Max number of snapshots %d" % N_max
    
    #print ''
    #for m in metadata:
    #    print m
    #print ''
    
    # now allocate the memory for the arrays
    K = len(metadata)
    
    T_k = numpy.ones(K,float)*options.temperature # inital temperatures are all equal
    beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
    
    #N_k = numpy.zeros([K], numpy.int32) # N_k[k] is the number of snapshots from umbrella simulation k
    # replaced by metadata[k]['n']
    
    #K_k = numpy.zeros([K], numpy.float64) # K_k[k] is the spring constant (in kJ/mol/deg**2) for umbrella simulation k
    # replaced by metadata[k]['k']
    
    #chi0_k = numpy.zeros([K], numpy.float64) # chi0_k[k] is the spring center location (in deg) for umbrella simulation k
    # replaced by metadata[k]['coord']
    
    #chi_kn = numpy.zeros([K,N_max], numpy.float64) # chi_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k
    data = numpy.zeros([K,N_max], numpy.float64) # the data!
    
    # TODO: combine the data somehow
    u_kn = numpy.zeros([K,N_max], numpy.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
    u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l
    
    g_k = numpy.zeros([K],numpy.float32) # statistical inefficiency?

    data_min = None
    data_max = None
    
    # Now loop through each datafile and extract the data
    for i, w in enumerate(metadata):
        print "Reading %s..." % w['path']
        f = open(w['path'], 'r')
        lines = f.readlines()
        f.close()
        
        clean_split_lines = [ filter(None, line.strip().expandtabs().split(' ')) for line in lines if not line.startswith('#') ]

        if different_temperatures:
            # if different temperatures are specified the metadata file, 
            # then we need the energies to compute the PMF, found in the third column
            for j,l in enumerate(clean_split_lines):
                data[i,j] = float(l[1]) # second column is the coordinate
                # third column will be the system's potential energy
                potential_energy = float(l[2])
                dchi = w['coord']-float(l[1])
                restraint_potential = 0.5*w['k']*(dchi**2)
                # TODO: given the coordinate and the restraining potential, calculate the umbrella restraint
                u_kn[i,j] = beta_k[i] * (potential_energy-restraint_potential) # reduced potential energy without umbrella restraint
        
            # Compute correlation times for potential energy and timeseries.
            # If the temperatures differ, use energies to determine samples; otherwise, use the cosine of chi
            g_k[i] = timeseries.statisticalInefficiency(u_kn[i,:], u_kn[i,:])
            indices = timeseries.subsampleCorrelatedData(u_kn[i,:])
        else:
            # no temperature column
            # FORMAT: f.write('0.0    %f\n' % deg)
            for j,l in enumerate(clean_split_lines):
                data[i,j] = float(l[1])
            dataset = numpy.cos(data[i,:w['n']])
            g_k[i] = timeseries.statisticalInefficiency(dataset,dataset)
            indices = timeseries.subsampleCorrelatedData(dataset)
            #g_k[i] = timeseries.statisticalInefficiency(numpy.cos(data[i,:]/(180.0/numpy.pi)),numpy.cos(data[i,:]/(180.0/numpy.pi)))
            #indices = timeseries.subsampleCorrelatedData(numpy.cos(data[i,:]/(180.0/numpy.pi)))

        # get min and max for data, used for binning ranges
        max = numpy.max(data[i,indices])
        min = numpy.min(data[i,indices])
        if data_max is None or data_max < max:
            data_max = max
        if data_min is None or data_min > min:
            data_min = min

        print "Correlation time for set %5d is %10.3f" % (i,g_k[i])

        # Subsample the data
        w['n'] = len(indices)
        u_kn[i,0:w['n']] = u_kn[i,indices]
        data[i,0:w['n']] = data[i,indices]

    print "Finished reading data files"
    
    # Set zero of u_kn -- this is arbitrary.
    u_kn -= u_kn.min()

    # Construct torsion bins
    print "Binning data..."
    print "Min coord: %f" % data_min
    print "Max coord: %f" % data_max
    
    delta = (data_max - data_min) / float(options.bins)
    print "Delta for binning %f" % delta
    # compute bin centers
    bin_center_i = numpy.zeros([options.bins], numpy.float64)
    for i in range(options.bins):
        bin_center_i[i] = data_min + delta/2 + delta * i
    #print "Bin centers:"
    #print bin_center_i
    
    # Bin data
    bin_kn = numpy.zeros([K,N_max], numpy.int32)
    # for each window
    for k in range(K):
        # for 0 to the number of snapshots in the window k
        for n in range(metadata[k]['n']):            
            # Compute bin assignment.
            bin_kn[k,n] = int((data[k,n] - data_min) / delta)
            for l in range(K):
                # Compute minimum-image torsion deviation from umbrella center l
                dchi = data[k,n] - metadata[l]['coord']
                # Compute energy of snapshot n from simulation k in umbrella potential l
                u_kln[k,l,n] = u_kn[k,n] + beta_k[k]*0.5*metadata[l]['k']*(dchi**2)
    
    for i in range(options.bins):
        print "bin %d" % i
        print numpy.sum(bin_kn==i)

    # Initialize MBAR.
    print "Running MBAR..."
    N_k = numpy.array([ w['n'] for w in metadata ], numpy.int32)
    print N_k
    print u_kln
    mbar = pymbar.MBAR(u_kln, N_k, verbose=True, method='self-consistent-iteration', initialize='BAR')
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson')

    # Compute PMF in unbiased potential (in units of kT).
    (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, options.bins)

    # Write out PMF
    print "PMF (in units of kT)"
    print "%8s %8s %8s" % ('bin', 'f', 'df')
    for i in range(options.bins):
        print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])
    
if __name__=='__main__':
    main()

