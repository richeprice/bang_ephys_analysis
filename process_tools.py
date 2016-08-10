#!/usr/bin/python
'''tools for processing raw data'''
import numpy as np
import scipy.cluster.vq as km

def rising_phases(trace,delta = .0002,max_delta = .05,min_run = 3):
    ''' find regions in trace with rising phase greater than delta
    returns: trace gradient and spike list of (start of rising phase,length of rising phase)'''
    tr_grad = np.gradient(trace)
    rising_ind = np.where(np.logical_and((tr_grad > delta),tr_grad < max_delta))[0]
    rising_diff = np.diff(rising_ind)
    #print len(rising_diff)
    is_rising = 0
    start = 0
    spike_list = []
    for i,v in enumerate(rising_diff):
        #print is_rising
        if v == 1:
            is_rising +=1
            if is_rising == min_run:
                start = rising_ind[i]
        else:
            is_rising = 0
            if start:
                if 4 < (rising_ind[i] - start) < 80:
                    spike_list.append([start,rising_ind[i]])
            start = 0
    return tr_grad,spike_list



def wiggle(positions,raw_data,win = 2,max_adjust = 20,x_direction = -1,y_direction = -1,threshold = .0008):
        ''' adjust spikes/stims toward peak '''
        if x_direction == -1:
            left = win
            right = 0
        elif x_direction == 1:
            left = 0
            right = win
        else:
            left = win/2
            right = win/2

        if y_direction == 1:
            extrema_func = np.max
            extrema_pos_func = np.argmax
        elif y_direction == -1:
            extrema_func = np.min
            extrema_pos_func = np.argmin

        adjusted_values = []
        for v in positions:
            new_v = v
            region = raw_data[v-left:v+right]

#            print "old pos = %s extrema = %s" %(raw_data[new_v],extrema_func(region))

            while (abs(raw_data[new_v] - extrema_func(region)) > threshold):
                #print "finding new_v"
                new_v = new_v - win + extrema_pos_func(region)
                region = raw_data[new_v - win: new_v + win]
                if len(region) == 0:
                    break

            if abs(new_v - v) > max_adjust:
                print 'adjusting too much!!'
                new_v = v
            adjusted_values.append(new_v)
            
        return(adjusted_values)


def compute_features_from_spike_mat(spike_mat):
    '''given a spike matrix where each row is a trace corresponding to spike, compute features and return a features matrix with columns as features'''

    maxes = []
    up_slopes = []
    down_slopes = []
    max_indices = []
    for trace in spike_mat:
        if len(trace) == 0:
            max_v = 1
            print "no trace"
        else:
            max_v = trace[10:35].max()
        if max_v > .15:
            print "too big"
            up_slopes.append(-10.)
            down_slopes.append(-10.)
            maxes.append(100)
            max_indices.append(100)
            continue
        grad = np.gradient(trace[10:30])
        max_i = trace[10:35].argmax()
        up_slope = (max_v - trace[10]) / max_i 
        down_slope = (trace[max_i + 20] - max_v) / 10.
        maxes.append(max_v)
        max_indices.append(grad[5])
        up_slopes.append(up_slope)
        down_slopes.append(down_slope)
    features_mat = np.array([np.array(max_indices),np.array(up_slopes),np.array(maxes)]).T  #[fwhm_ar.np.array(up_slopes),np.array(down_slopes),np.array(maxes)]).T
    return(features_mat)

def cluster_neurons(features_mat,positions,k = 6):
    '''given a features matrix and positions list, return a dict of lists with each key being a neuron code (0...k-1) and values be a list of global spike positions for each neuron'''
    neurons_dic = {}
    obs = km.whiten(features_mat)
    codebook,d = km.kmeans(obs,k)
    codes,dist = km.vq(obs,codebook)
    for i,code in enumerate(codes):
        try:
            neurons_dic[code].append(positions[i])
        except:
            neurons_dic[code] = [positions[i]]
    for code in neurons_dic:
        neurons_dic[code] = np.array(neurons_dic[code])
    return(neurons_dic)
