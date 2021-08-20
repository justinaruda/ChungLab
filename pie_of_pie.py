import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np
import pandas as pd

def plot_pies(name,table_df,loci_type):
    destabilized_pcts = dict()
    stabilized_pcts = dict()
    no_change_pcts = dict()
    
    destabilized = table_df[table_df['Raw_FC']<=0.95]
    destabilized_pct = len(destabilized.index.values)/len(table_df.index.values)
    stabilized = table_df[table_df['Raw_FC']>=1.05]
    stabilized_pct = len(stabilized.index.values)/len(table_df.index.values)
    no_change = table_df[(table_df['Raw_FC']>0.95) & (table_df['Raw_FC']<1.05)]
    no_change_pct = len(no_change.index.values)/len(table_df.index.values)
    
    destabilized_ir_clusters = destabilized[destabilized['Repeat_Context']=='ir_cluster']
    destabilized_pcts['ir'] = len(destabilized_ir_clusters.index.values)/len(destabilized.index.values)
    destabilized_clusters = destabilized[destabilized['Repeat_Context']=='cluster']
    destabilized_pcts['cluster'] = len(destabilized_clusters.index.values)/len(destabilized.index.values)
    destabilized_individual = destabilized[destabilized['Repeat_Context']=='individual']
    destabilized_pcts['individual'] = len(destabilized_individual.index.values)/len(destabilized.index.values)
    destabilized_none = destabilized[destabilized['Repeat_Context']=='None']
    destabilized_pcts['none'] = len(destabilized_none.index.values)/len(destabilized.index.values)
    
    stabilized_ir_clusters = stabilized[stabilized['Repeat_Context']=='ir_cluster']
    stabilized_pcts['ir'] = len(stabilized_ir_clusters.index.values)/len(stabilized.index.values)
    stabilized_clusters = stabilized[stabilized['Repeat_Context']=='cluster']
    stabilized_pcts['cluster'] = len(stabilized_clusters.index.values)/len(stabilized.index.values)
    stabilized_individual = stabilized[stabilized['Repeat_Context']=='individual']
    stabilized_pcts['individual'] = len(stabilized_individual.index.values)/len(stabilized.index.values)
    stabilized_none = stabilized[stabilized['Repeat_Context']=='None']
    stabilized_pcts['none'] = len(stabilized_none.index.values)/len(stabilized.index.values)
    
    no_change_ir_clusters = no_change[no_change['Repeat_Context']=='ir_cluster']
    no_change_pcts['ir'] = len(no_change_ir_clusters.index.values)/len(no_change.index.values)
    no_change_clusters = no_change[no_change['Repeat_Context']=='cluster']
    no_change_pcts['cluster'] = len(no_change_clusters.index.values)/len(no_change.index.values)
    no_change_individual = no_change[no_change['Repeat_Context']=='individual']
    no_change_pcts['individual'] = len(no_change_individual.index.values)/len(no_change.index.values)
    no_change_none = no_change[no_change['Repeat_Context']=='None']
    no_change_pcts['none'] = len(no_change_none.index.values)/len(no_change.index.values)
    
    
    # make figure and assign axis objects
    
    fig = plt.figure(figsize=(16, 10))
    ax1 = fig.add_subplot(232)
    ax2 = fig.add_subplot(233)
    ax3 = fig.add_subplot(235)
    ax4 = fig.add_subplot(231)
    fig.subplots_adjust(wspace=-0.3,hspace=-0.3)
    
    theme = plt.get_cmap('Pastel1')
    
    ax1.set_prop_cycle("color", [theme(1. * i / 3)
                                 for i in range(3)])
    
    theme = plt.get_cmap('viridis')
    ax2.set_prop_cycle("color", [theme(1. * i / 4)
                                 for i in range(4)])
    ax3.set_prop_cycle("color", [theme(1. * i / 4)
                                 for i in range(4)])
    ax4.set_prop_cycle("color", [theme(1. * i / 4)
                                 for i in range(4)])
    
    # large pie chart parameters
    ratios = [destabilized_pct,no_change_pct,stabilized_pct]
    labels = ['Destabilized','No Change','Stabilized']
    explode = [0.1, 0.1, 0.1]
    # rotate so that first wedge is split by the x-axis
    angle = (-180 * ratios[0]) - 40
    ax1.pie(ratios, autopct='%1.1f%%', startangle=angle, explode=explode)
    
    # small pie chart parameters
    destabilized_ratios = [destabilized_pcts['ir'],destabilized_pcts['cluster'],destabilized_pcts['individual'],destabilized_pcts['none']]
    destabilized_labels = [100*x for x in destabilized_ratios]
    stabilized_ratios = [stabilized_pcts['ir'],stabilized_pcts['cluster'],stabilized_pcts['individual'],stabilized_pcts['none']]
    stabilized_labels = [100*x for x in stabilized_ratios]
    no_change_ratios = [no_change_pcts['ir'],no_change_pcts['cluster'],no_change_pcts['individual'],no_change_pcts['none']]
    no_change_labels = [100*x for x in no_change_ratios]
    
    labels = ['IR Alu Cluster', 'Alu Cluster', 'Individual Alu', 'Non-Alu']
    width = .2
    
    ax2.pie(destabilized_ratios, startangle=angle, radius=0.5, textprops={'size': 'smaller'})
    labels2 = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(labels, destabilized_labels)]
    
    ax3.pie(stabilized_ratios, startangle=angle, radius=0.5, textprops={'size': 'smaller'})
    labels3 = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(labels, stabilized_labels)]
    
    ax4.pie(no_change_ratios, startangle=angle, radius=0.5, textprops={'size': 'smaller'})
    labels4 = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(labels, no_change_labels)]
    
    ax1.set_title('Predicted $\Delta_{{Editing}}$ Stability of 1000nts Flanking {0} {1}'.format(name.replace('_',' '),loci_type),size='x-large')
    ax2.set_title('Destabilized',y=0.01)
    ax3.set_title('Stabilized',y=0.01)
    ax4.set_title('No Change',y=0.01)
    
    ax2.legend(labels=labels2,loc='lower center',bbox_to_anchor=(0.5,0.065))
    ax3.legend(labels=labels3,loc='lower center',bbox_to_anchor=(0.5,0.065))
    ax4.legend(labels=labels4,loc='lower center',bbox_to_anchor=(0.5,0.065))
    
    # use ConnectionPatch to draw lines between the two plots
    # get the wedge data
    
    for ax,patch in zip([ax2,ax4,ax3],ax1.patches):
        
        theta1, theta2 = patch.theta1, patch.theta2
        center, r = patch.center, patch.r
        x = r * np.cos(np.pi / 180 * theta2) + center[0]
        y = np.sin(np.pi / 180 * theta2) + center[1]
        con = ConnectionPatch(xyA=(- width / 2, 0), xyB=(x, y),coordsA="data", coordsB="data", axesA=ax, axesB=ax1)
        con.set_color([0, 0, 0])
        con.set_zorder(-1)
        con.set_linewidth(2)
        ax2.add_artist(con)
        
    
    # destabilized_theta1, destabilized_theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
    # center, r = ax1.patches[0].center, ax1.patches[0].r
    
    # stabilized_theta1, destabilized_theta2 = ax1.patches[1].theta1, ax1.patches[1].theta2
    # center, r = ax1.patches[1].center, ax1.patches[1].r
    
    # no_change_theta1, no_change_theta2 = ax1.patches[2].theta1, ax1.patches[2].theta2
    # center, r = ax1.patches[2].center, ax1.patches[2].r
    
    # draw top connecting line
    # x = r * np.cos(np.pi / 180 * theta2) + center[0]
    # y = np.sin(np.pi / 180 * theta2) + center[1]
    # con = ConnectionPatch(xyA=(- width / 2, .5), xyB=(x, y),coordsA="data", coordsB="data", axesA=ax2, axesB=ax1)
    # con.set_color([0, 0, 0])
    # con.set_linewidth(2)
    # ax2.add_artist(con)
    
    # # draw bottom connecting line
    # x = r * np.cos(np.pi / 180 * theta1) + center[0]
    # y = np.sin(np.pi / 180 * theta1) + center[1]
    # con = ConnectionPatch(xyA=(- width / 2, -.5), xyB=(x, y), coordsA="data",coordsB="data", axesA=ax2, axesB=ax1)
    # con.set_color([0, 0, 0])
    # ax2.add_artist(con)
    # con.set_linewidth(2)
#    plt.savefig(r'C:/Users/justi/Desktop/RNAfold/figures/{0}_{1}.png'.format(name,loci_type))
    plt.show()
    
    fig = plt.figure(figsize=(16, 10))
    ax1 = fig.add_subplot(111)
    theme = plt.get_cmap('viridis')
    ax1.set_prop_cycle("color", [theme(1. * i / 4)
                                 for i in range(4)])
    
    
    ir_alus = table_df[table_df['Repeat_Context']=='ir_cluster']
    ir_alus_pct = len(ir_alus.index.values)/len(table_df.index.values)
    cluster = table_df[table_df['Repeat_Context']=='cluster']
    cluster_pct = len(cluster.index.values)/len(table_df.index.values)
    individual = table_df[table_df['Repeat_Context']=='individual']
    individual_pct = len(individual.index.values)/len(table_df.index.values)
    none = table_df[table_df['Repeat_Context']=='None']
    none_pct = len(none.index.values)/len(table_df.index.values)
    
    ax1.set_title('Alu Containing Sequence Breakdown {0} {1}'.format(name.replace('_',' '),loci_type),size='xx-large')
    pcts = [ir_alus_pct,cluster_pct,individual_pct,none_pct]
    pct_labels = [100*x for x in pcts]
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(['IR Alu Cluster', 'Alu Cluster', 'Individual Alu', 'Non-Alu'],pct_labels)]
    angle = (-180 * ratios[0]) - 40
    ax1.pie(pcts,startangle=angle)
    ax1.legend(labels=labels,loc='lower center',bbox_to_anchor=(0.5,-0.08),fontsize='x-large')
    plt.savefig(r'C:/Users/justi/Desktop/RNAfold/figures/{0}_{1}_Alu_Breakdown.png'.format(name,loci_type))
    plt.show()
    
    
    
    
    
#table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/CITS/A1KO_Mock_Down.pool.tag.uniq.del.CITS.s30.singleton.1001nt.unedited_vs_A1KO_Mock_Down.pool.tag.uniq.del.CITS.s30.singleton.1001nt.edited_compared.ir.out')

for table in ['A1KO_Mock_Down','A1KO_IFN_Down','A1KO_Mock_Up','A1KO_IFN_Up']:
    CITS_table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/CITS/annotated/{0}.pool.tag.uniq.del.CITS.s30.singleton.1001nt.unedited_vs_{0}.pool.tag.uniq.del.CITS.s30.singleton.1001nt.edited_compared.ir.out'.format(table))
    plot_pies(table,CITS_table_df,'CITS')
#    peaks_table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/peaks/annotated/{0}.pool.tag.uniq.peak.sig.1001nt.unedited_vs_{0}.pool.tag.uniq.peak.sig.1001nt.edited_compared.ir.out'.format(table))    
#    plot_pies(table,peaks_table_df,'Peaks')

#table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/CITS/random_intervals.1001nt.unedited_vs_random_intervals.1001nt.edited_compared.ir.out')

# table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/CITS/old/random_intervals.1001nt.unedited_vs_random_intervals.1001nt.edited_compared.ir.out')
# plot_pies('Random_old',table_df,'Loci')

table_df = pd.read_table(r'C:/Users/justi/Desktop/RNAfold/CITS/annotated/random_intervals.1001nt.unedited_vs_random_intervals.1001nt.edited_compared.ir.out')
plot_pies('Random',table_df,'Loci')