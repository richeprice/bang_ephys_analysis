#!/usr/bin/python

''' collection of tools for analysis of spike train times'''
import pylab as py
import seaborn as sns
import numpy as np
def compute_firing_rate():
    pass
exp_decay = lambda x, A1, t1, y0: A1 * np.exp(x * t1) + y0
rgb = [['255', '127', '14'],
       ['214', '39', '40'],
       ['148', '103', '189'],
       ['44', '160', '44'],
       ['31', '119', '180'],
       ['227', '119', '194'],
       ['188', '189', '34'],
       ['140', '86', '75'],
       ['127', '127', '127'],
       ['23', '190', '207']]

sns.set_style({'legend.frameon':True})

def compute_coherence(st1,st2,freq=np.linspace(0,200,100000)):
    ft1 = np.array([np.sum(np.exp(-2j*k*np.pi*st1)) for k in freq])
    ft2 = np.array([np.sum(np.exp(-2j*k*np.pi*st2)) for k in freq])
    st1_squ = ft1**2
    st2_squ = ft2**2
    

def pull_ecs_spikes(st,ecs,down = 22050 * 50):
    '''given a spiketrain and array/list of ecs times, return a list of spike train arrays for each ecs, spike times relative to ecs position'''
    st_list = []
    for pos in ecs:
        start = st.searchsorted(pos)
        end = st.searchsorted(pos+down) - 1
        st_list.append(st[start:end] - pos)
    return(st_list)

def window_data(data,window_size,func):
    pass

def cross_correlation(a,b,up,down):
    ''' given arrays a and b, return crosscorrelation of a in b'''
    window = np.zeros(up+down)
    for pos_a in a:
        start = b.searchsorted(pos_a - up)
        end = b.searchsorted(pos_a + down)
        for pos_b in b[start:end]:
            window[pos_b - pos_a - up]+=1
    return (window)

def latency_by_last_spike(a,b):
    ''' take spikes in a and find latency from last b and last a spike '''
    li = []
    for i,pos_a in enumerate(a[1:]):
        last_b = b[b.searchsorted(pos_a) - 1]
        last_a = a[i]
        li.append([pos_a-last_b,pos_a-last_a])
    return np.array(li)

def latency_by_prev_GF(a,b):
    ''' take spikes in a and find latency from last b and last a spike '''
    li = []
    for i,pos_a in enumerate(a[1:]):
        last_b = b[b.searchsorted(pos_a) - 1]
        last_a = a[i]
        li.append([pos_a-last_b,last_b - b[b.searchsorted(pos_a)-2]])
    return np.array(li)


def firing_rate_at_positions(positions,st,forward = 22050,back = 22050):
    '''for position in positions, determine firing rate from spike times(st) '''
    fr = np.zeros(len(positions))
    for i,position in enumerate(positions):
        fr[i] = float(st.searchsorted(position + forward) - st.searchsorted(position - back))/(back+forward)
    return fr




def gen_gaussian(length,sigma = .5,boundaries = (-2,2)):
    ''' return array of weights for adjustment'''
    c = np.sqrt(sigma)
    b = 0
    a = 1/(c*np.sqrt(2*3.141592654))
    ar = np.zeros(length)
    for i,v in enumerate(np.linspace(boundaries[0],boundaries[1],length)):
        ar[i] = gaussian(v,a,b,c)
    return ar

def gen_exp(length,tau):
    return np.array([exp_decay(x,1,-.5,0) for x in np.linspace(0,10,length)])



def stack_spikes(spike_times,length = 0,win_size = 2*22050):
    #    stack_ar = gen_gaussian(win_size)
    stack_ar = gen_exp(win_size,-.5)
    
    new_ar = np.zeros(length)
    for spike in spike_times:
        try:
             new_ar[spike:spike+win_size] += stack_ar
        except:
            continue
    return(new_ar)


def resize_ar(ar,scaling,window = None):
    '''given an array of continous data , return a resized array scaled by scaling '''
    if not(window):
        window = scaling*2
    new_ar = np.zeros((len(ar)/scaling))
    for i,v in enumerate(np.linspace(window/2,len(ar) - window/2,len(new_ar))):
        new_ar[i] = ar[int(v)-window/2:int(v)+window/2].mean()
    return new_ar
                       
    


############# epoch/meta tools

def remove_stimulus_responses(m_train,stims,window = 300,return_responses = False):
    '''given a m spike train, remove spikes within window of stim, return new m_train,  return responses list if return_respones = True'''
    removal_indices = []
    for train in m_train:
        removal_indices.append(range(len(train)))
    responses = {}
    for stim in stims:
        responses[stim] = []
        for i,train in enumerate(m_train):
            spike_i = train.searchsorted(stim)
            if spike_i == len(train):
                responses[stim].append(None)
                continue
            latency = train[spike_i] - stim
            if latency < window:
                try:
                    removal_indices[i].remove(spike_i)
                except(ValueError):
                    print "ch %s has removed index %s" %(i,spike_i)
                if return_responses:
                    responses[stim].append(latency)
            elif return_responses:
                responses[stim].append(None)
    new_m_train = [m_train[i][removal_indices[i]] for i in xrange(len(m_train))]
    return(new_m_train,responses)






def list_of_st_to_fr(lst,padding = 22050):
    lfr = []
    for st in lst:
        length = st[0][-1] - st[0][0] + padding
        adjusted_st = np.array([pos+padding/2 for pos in st[0]])
        lfr.append(stack_spikes(adjusted_st,length,win_size = 22050))

    return(lfr)

def list_of_m_st_to_fr(lst,padding = 22050):
    li_mfr = []
    for m_st in lst:
        m_fr = []
        for st in m_st:
            if len(st) == 0:
                m_fr.append([])
                continue
            length = st[-1] - st[0] + padding
            adjusted_st = np.array([pos+padding/2 for pos in st])
            m_fr.append(stack_spikes(adjusted_st,length,win_size = 22050))
        li_mfr.append(m_fr)
    return(li_mfr)


def unwrap_stims_arr(stims_arr):
    li = []
    for ar in stims_arr:
        li.extend(ar)
    ar = np.array(li)
    ar.sort()
    return(ar)
        
def sort_ecs(stims_dic):
    return(np.array(stims_dic['ecs_n']).argsort())

def generate_stims_dic(stims_arr):
    stim_dict = {'ecs_n':[]}
    for stims in stims_arr:
        if len(stims) == 1:
            key = 's'
        elif len(stims) ==2:
            key = 'pp'
        elif len(stims) == 10:
            key = 'ff'
        elif len(stims) > 15:
            stim_dict['ecs_n'].append(len(stims))
            key = 'ecs'
        else:
            key = 'unk'
        try:
            stim_dict[key].append(stims[-1])
        except:
            stim_dict[key] = [stims[-1]]
    return(stim_dict)
            

def generate_header_dic(header):
    ''' given a header "geno_bottle_sex_age_time" return a dictionary'''
    keys = ['geno','bottle','sex','age','time']
    dic = dict(zip(keys,['']*len(keys)))
    split_head = header.split("_")
    if len(split_head) != 5:
        print "wierd head"
        return(dic)
    dic = dict(zip(keys,split_head))
    return dic

def pull_spikes_by_epoch(st,start,length):
    ''' return spikes from spike train (st) within epoch, return times relative to epoch start '''
    ar = st[st.searchsorted(start):st.searchsorted(start+length)] - start
    return ar



def collect_sub_spike_train(st,min_dist,min_spikes = 3,as_array = True):
    ''' given a 1d spike train, return a list of spike trains of min(min_spikes) spikes and greater than min_dist(in samples) between trains'''
    padded_st = [-1*min_dist]   ## changed from 0 to -1*min_dist
    padded_st.extend(st)
    padded_st.append(st[-1] + min_dist + 1)
    padded_st_ar = np.array(padded_st)
    end_spikes = np.where(np.diff(padded_st[1:]) > min_dist)[0] # indices in st where the i and i+1 are greater than min_dist,  ith position is end of a train
    start_spikes = np.where(np.diff(padded_st[:-1]) > min_dist)[0]
    st_li = []
    positions_li = []
    for s in start_spikes:
        start_time = st[s]
        e = end_spikes[end_spikes.searchsorted(s)]
        end_time = st[e]
        if e-s < min_spikes:
            continue
        st_li.append(sub_mspike_train([st],start_time,end_time,as_array = as_array))
        positions_li.append([start_time,end_time])
    return(st_li,positions_li)


def find_bursts(st,min_dist,min_spikes,as_array = True):
    
    gaps = np.where(np.diff(st) > min_dist)[0]
    start_ends = []
    start = 0
    
    for gap in gaps:
        end = gap + 1
        start_ends.append((start,end))
        start = gap + 1
    start_ends.append((start,len(st)))

    st_li = []
    pos_li = []
    for s,e in start_ends:
        if e-s > min_spikes:
            st_li.append(st[s:e])
#            pos_li.append((st[s],st[e]))
    return(st_li)
            



def collect_sub_m_spike_train(m_st,min_dist,min_spikes = 3,as_array = True):
    ''' given a m spike train, return a list of spike trains of min(min_spikes) spikes and greater than min_dist(in samples) between trains'''
    st = combine_m_spike_train(m_st)
    padded_st = [-1*min_dist]
    padded_st.extend(st)
    padded_st.append(st[-1] + min_dist + 1)
    padded_st_ar = np.array(padded_st)
    end_spikes = np.where(np.diff(padded_st[1:]) > min_dist)[0] # indices in st where the i and i+1 are greater than min_dist,  ith position is end of a train
    start_spikes = np.where(np.diff(padded_st[:-1]) > min_dist)[0]
    st_li = []
    positions_li = []
    for s in start_spikes:
        start_time = st[s]
        e = end_spikes[end_spikes.searchsorted(s)]
        end_time = st[e]
        if e-s < min_spikes:
            continue
        st_li.append(sub_mspike_train(m_st,start_time,end_time,as_array = as_array))
        positions_li.append([start_time,end_time])
    return(st_li,positions_li)

def combine_m_spike_train(m_spikes):
    ''' given multi spike train, return single spike train with combined spikes'''
    all_spikes = []
    for st in m_spikes:
        all_spikes.extend(st)
    st = np.array(all_spikes)
    st.sort()
    return(st)

def m_train_to_obs(m_train):
    return([[[x for x in ar] for ar in m_train]])



def sub_mspike_train(m_train,start,end,relative = True,as_array = True):
    ''' given a multichannel spike train [st0,st1,...stn], return a sub m_train bounded by start,end'''
    sub_mtrain = []
    for train in m_train:
        sub_train = train[train.searchsorted(start):train.searchsorted(end)]
        if relative:
            sub_train = [pos - start for pos in sub_train]
        sub_mtrain.append(sub_train)
    if as_array:
        return([np.array(train) for train in sub_mtrain])
    else:
        return(sub_mtrain)



def pad_li_li_arrays(li_li_arrays,length = 0):
    ''' given a list of list of arrays, find the length of the longest array (max_length) and pad all arrays to length max_length
    returns li of li of padded arrays and max_length'''
    return_li = []
    if not(length):
        length = np.max([len(ar) for ar in np.ravel(li_li_arrays)])
    num_chs = len(li_li_arrays[0])
    for trial  in li_li_arrays:
        trial_li = []
        for ar in trial:
            data = np.zeros(length)
            data[:len(ar)] = ar
            trial_li.append(data)
        return_li.append(trial_li)
    return(return_li,length)
        
        
                         




def concatenate_li_li_ar(m_spikes_li):
    ''' given a list of m_spikes, return single m_spike with concatenated spiketrains'''
    ch_dic = {}
    ch_li = []
    for chs in m_spikes_li:
        for i,ch in enumerate(chs):
            try:
                ch_dic[i].append(ch)
            except:
                ch_dic[i] = [ch]
    for i in xrange(len(ch_dic.keys())):
        ch_li.append(np.concatenate(ch_dic[i]))
    return ch_li
            

def generate_ecs_dic(stims_arr,sizes,ecs_li):
    ''' given stims_arr (array of arrays with variable length),sizes = len of concatenated stims_arr,ecs_li = list of times in stims_arr for last stim of ecs
    returns a dictionary keyed by ecs_li times, values = array of len(num_stims) with values from sizes'''
    ecs_dic = {}
    sizes_pos = 0
    for stims in stims_arr:
        if stims[-1] in ecs_li:
            ecs_dic[stims[-1]] = sizes[sizes_pos:sizes_pos + len(stims)]
        sizes_pos += len(stims)
    return(ecs_dic)

def get_post_ecs_stims(stims_dic,downstream = 60):
    ecs_li = stims_dic['ecs']
    ss = np.array(stims_dic['s'])
    stim_coords_li = []
    for ecs in ecs_li:
        stim_coords_li.append(ss[ss.searchsorted(ecs):ss.searchsorted(ecs + 60 * 22050)])
    return stim_coords_li




### plotting tools




def normalize_list_of_chs(li_of_chs):
    '''rescale values in each ch array to 0 to 1'''
    normalized_li = []
    chs_all = concatenate_li_li_ar(li_of_chs)
    chs_max = [float(ar.max()) for ar in chs_all]
    for chs in li_of_chs:
        new_chs = []
        for i,ch in enumerate(chs):
            new_chs.append(ch/chs_max[i])
        normalized_li.append(new_chs)
    return(normalized_li)




def plot_list_of_chs(li_ch,ymax = .8,xtime = None,header_li = None):
    '''make a subplot with traces from each channel '''
    backgrounds = ['whitesmoke','gainsboro']
    num_subplots = len(li_ch)
    subplots_li = []
    fig,subplots = py.subplots(num_subplots,1,sharex = True,sharey = True)
    if num_subplots == 1:
        subplots_li = [subplots]
    else:
        subplots_li.extend(subplots)
    for i,(li,ax) in enumerate(zip(li_ch,subplots_li),1):
        ax.set_axis_bgcolor(backgrounds[i%2])
        ax.set_ylim(0,ymax)
        ax.set_yticklabels([])
        end = 0
        for ch in li:
            end = len(ch)
            if xtime:

                xs = np.linspace(0,xtime,end)
                ax.plot(xs,ch,linewidth = 3)
                end = xtime
            else:
                ax.plot(ch,linewidth = 3)
        if header_li:
            ax.text(end - 5,ymax - .2,header_li[i-1],fontsize = 8)
#            ax.text(end - 100,ymax - 2,header_li[i-1],fontsize = 8)
        ax.text(end-10,2,str(i-1))
        ax.set_xticks(np.arange(0,xtime,5))
        ax.set_yticks([])
#        ax.set_yscale('log')
        ax.grid(True)


    fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.95)
    
    return(fig)

def plot_list_of_chs_with_stimulus(li_ch,ymax = 5,xtime = None,header_li = None):
    '''make a subplot with traces from each channel '''
    backgrounds = ['whitesmoke','gainsboro']
    num_subplots = len(li_ch)
    subplots_li = []
    fig,subplots = py.subplots(num_subplots,1,sharex = True)
    if num_subplots == 1:
        subplots_li = [subplots]
    else:
        subplots_li.extend(subplots)
    for i,(li,ax) in enumerate(zip(li_ch,subplots_li),1):
        ax.set_axis_bgcolor(backgrounds[i%2])
        ax.set_ylim(0,ymax)
        ax.set_yticklabels([])
        end = 0
        for j,ch in enumerate(li):
            end = len(ch)
            if xtime:
                end = xtime
                xs = np.linspace(0,xtime,end)
                ax.plot(xs,ch,linewidth = 2,color = [float(c)/255 for c in rgb[j]])
            else:
                ax.plot(ch,color = [float(c)/255 for c in rgb[j]],linewidth = 2)
        if header_li:
            print 'header'
#            ax.text(end - 5,ymax - .2,header_li[i-1],fontsize = 8)
            ax.text(0,0,header_li[i-1],fontsize = 8)
        ax.text(end-10,2,str(i-1))
        ax.set_xticks(np.arange(0,xtime,5))
        ax.set_yticks([])
        ax.grid(True)

    ax.set_xlabel('Time (sec)',fontsize = 18)
    fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.95)
    
    return(fig)


def plot_traces_from_channels(epoch_traces,time,si,ch_list = [0,3,1,2],spike_trains = 0,firing_rates = 0,step = 5,alternate_bg = True,offset = 0,labels=['DVM-r','DLM-r','DLM-l','DVM-l']):
    palette = sns.palettes.SEABORN_PALETTES['colorblind']
    yloc = [0.92,.90,.86,.83]
    pal = [palette[x] for x in [2,0,1,3,4,5]]
    num_channels = len(ch_list)
    fig,axes = py.subplots(len(si)*len(ch_list),1,sharex=True,sharey=True)
    for i,x in enumerate(si):
        [ax.plot(time,epoch_traces[x][ch][:len(time)],color = pal[ch],label = labels[ch]) for j,(ax,ch) in enumerate(zip(axes[num_channels*i:i*num_channels+num_channels],ch_list))]
        if spike_trains:
            print i,x
            [ax.broken_barh([(float(li[0])/22050.,float(li[-1]-li[0])/22050.) for li in spike_trains[x][ch]],[-10,100],color = pal[ch],alpha = 0.2) for ax,ch in zip(axes[num_channels*i:i*num_channels+num_channels],ch_list)]
        if firing_rates:
            fr_time = np.linspace(0,time[-1],len(firing_rates[0][0]))
#            [ax.plot(fr_time,firing_rates[x][ch],color = pal[ch],alpha = 0.6,linewidth = 4) for ax,ch in zip(axes[num_channels*i:i*num_channels+num_channels],ch_list)]
            [ax.plot(fr_time,firing_rates[x][ch],color = 'black',alpha = .6,linewidth = 2) for ax,ch in zip(axes[num_channels*i:i*num_channels+num_channels],ch_list)]
        

        if alternate_bg:  ## change bg colors for each trial
            if i%2 == 1:
                print "changing axis bg lightgrey"
                [ax.set_axis_bgcolor('gainsboro') for ax in axes[num_channels*i:i*num_channels+num_channels]]

#        [ax.set_ylabel(i+1+offset,rotation = 'horizontal') for ax in axes[num_channels*i:i*num_channels+num_channels]]
        axes[num_channels*i].set_ylabel((i+1+offset),rotation = 'horizontal',size = 18)
                

    ax = axes[-1]
    ax.set_ylim(-0.01,0.1)
    ax.set_yticks([])
    ax.set_xlabel('time (sec)',fontsize = 16)
    ax.set_xticks(np.arange(0,time[-1]+step,step))
    ax.set_xticklabels(np.arange(0,len(time)+step,step),fontsize = 14)
    if labels:
        lines = [py.mpl.lines.Line2D(np.arange(10),np.arange(10),lw = 15,color = pal[ch]) for ch in ch_list]
        fig.legend(lines,[labels[ch] for ch in ch_list],loc = [.88,yloc[num_channels-1]],fancybox=True,shadow=True,prop={'size':20})
    

    fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.97,left = 0.03,right = 0.98)
    return fig








def plot_traces_from_channels_adjacent(epoch_traces,time,si,ch_list = [0,3,1,2],spike_trains = 0,firing_rates = 0,step = 5,alternate_bg = True,offset = 0,labels=['DVM-r','DLM-r','DLM-l','DVM-l']):
    ''' given list of epochs with list of arrays, time array, order list, channel list, plot each channel adjacent '''
    palette = sns.palettes.SEABORN_PALETTES['colorblind']
    yloc = [0.92,.90,.86,.83]
    pal = [palette[x] for x in [2,0,1,3,4,5]]
    num_channels = len(ch_list)
    fig,axes = py.subplots(len(si),len(ch_list),sharex=True,sharey=True)
    for i,x in enumerate(si):
        [ax.plot(time,epoch_traces[x][ch][:len(time)],color = pal[ch],label = labels[ch]) for j,(ax,ch) in enumerate(zip(axes[i],ch_list))]
        if spike_trains:
            print i,x
            [axes[i][0].broken_barh([(float(li[0])/22050.,float(li[-1]-li[0])/22050.) for li in spike_trains[x][ch]],[-10,100],color = pal[ch],alpha = 1) for ax,ch in zip(axes[i],ch_list)]
            #[ax.broken_barh([(float(li[0])/22050.,float(li[-1]-li[0])/22050.) for li in spike_trains[x][ch]],[-10,100],color = pal[ch],alpha = 0.2) for ax,ch in zip(axes[i],ch_list)]
        if firing_rates:
            fr_time = np.linspace(0,time[-1],len(firing_rates[0][0]))
#            [ax.plot(fr_time,firing_rates[x][ch],color = pal[ch],alpha = 0.6,linewidth = 4) for ax,ch in zip(axes[num_channels*i:i*num_channels+num_channels],ch_list)]
            [ax.plot(fr_time,firing_rates[x][ch],color = 'black',alpha = .6,linewidth = 2) for ax,ch in zip(axes[i],ch_list)]
        

        if alternate_bg:  ## change bg colors for each trial
            if i%2 == 1:
                print "changing axis bg lightgrey"
                [ax.set_axis_bgcolor('gainsboro') for ax in axes[i]]

#        [ax.set_ylabel(i+1+offset,rotation = 'horizontal') for ax in axes[num_channels*i:i*num_channels+num_channels]]
        axes[i][0].set_ylabel((i+1+offset),rotation = 'horizontal',size = 18)
                

    ax = axes[0][-1]
    ax.set_ylim(-0.01,0.1)
    ax.set_yticks([])
    ax.set_xlabel('time (sec)',fontsize = 16)
    ax.set_xticks(np.arange(0,time[-1]+step,step))
    ax.set_xticklabels(np.arange(0,len(time)+step,step),fontsize = 14)
    if labels:
        lines = [py.mpl.lines.Line2D(np.arange(10),np.arange(10),lw = 15,color = pal[ch]) for ch in ch_list]
        fig.legend(lines,[labels[ch] for ch in ch_list],loc = [.88,yloc[num_channels-1]],fancybox=True,shadow=True,prop={'size':20})
    

    fig.subplots_adjust(hspace = 0,bottom = 0.05,top = 0.97,left = 0.03,right = 0.98)
    return fig




def plot_firing_rates_by_channel(epochs_li,si,ch_list = [0,3,1,2],time = 60):
    ''' takes list of epochs and sorting indices and produces a figure'''
    t_fr = np.linspace(0,time,len(epochs_li[0][0]))
    fig,axes = py.subplots(4,1,sharex = True,sharey = True)
    [[ax.plot(t_fr,epochs_li[x][j],color = cmap(float(i)/len(epochs_li))) for i,x in enumerate(si)] for j,ax in zip(ch_list,axes)]
    axes[-1].set_xlabel('time (sec)')
    return fig
