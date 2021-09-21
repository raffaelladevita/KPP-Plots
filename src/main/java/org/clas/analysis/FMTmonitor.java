/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.utils.groups.IndexedTable;


/**
 *
 * @author devita
 */
public class FMTmonitor extends AnalysisMonitor {

    double rfPeriod = 4.008;
    public static final int FVT_Nlayers = 6;
    public static int FVT_Nstrips = 1024;// Number of strips
    public static int FVT_Halfstrips = 320; // In the middle of the FMT, 320 strips are split in two. 
    public static double FVT_Pitch = 0.0525; //strip width
    public static double FVT_Beamhole = 4.2575;//Radius of the hole in the center for the beam.

    double[] FVT_stripsXlocref;
    double[] FVT_stripsYlocref;
    double[] FVT_stripslength; //Give the strip length
    double[][][] FVT_stripsX; //Give the  end-points x-coordinates of the strip segment rotated in the correct frame for the layer
    double[][][] FVT_stripsY; //Give the  end-points y-coordinates of the strip segment

    double[] FVT_Zlayer = {293.197, 305.097, 316.997, 330.897, 342.797, 354.697}; //Give z-coordinate of the layer
//    double[] FVT_Alpha = {19, 79, 139, -161, -101, -41}; //Give the rotation angle to apply
//    double[] FVT_Alpha = {157.5, 217.5, 277.5, 157.5, 217.5, 277.5}; //Give the rotation angle to apply
    double[] FVT_Alpha = {67.5, 127.5, 187.5, 67.5, 127.5, 187.5}; //Give the rotation angle to apply
    double[] MY_Alpha = {67.5, 127.5, 187.5, 67.5, 127.5, 187.5}; //Give the rotation angle to apply

    int ntraj=0;
    int nfmt=0;
    
    public FMTmonitor(String name, ConstantsManager ccdb) {
        super(name,ccdb);
        this.setAnalysisTabNames("Monte Carlo","Truth", "Vertex", "FMT hits", "FMT phi", "Clusters", "Cluster phi", "Cluster phi vs strip", "Cluster strip vs strip", "Cluster energy");
        this.init(false);
        this.LoadGeometry();
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary", "summary", 6, 1, 7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1, 1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);

        // mc comparison
        H1F hi_dp_pos = new H1F("hi_dp_pos", "hi_dp_pos", 100, -0.1, 0.1); 
        hi_dp_pos.setTitleX("#Delta P/P");
        hi_dp_pos.setTitleY("Counts");
        hi_dp_pos.setTitle("Positive Tracks");
        F1D f1_dp_pos = new F1D("f1_dp_pos","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dp_pos.setParameter(0, 0);
        f1_dp_pos.setParameter(1, 0);
        f1_dp_pos.setParameter(2, 1.0);
        f1_dp_pos.setLineWidth(2);
        f1_dp_pos.setLineColor(2);
        f1_dp_pos.setOptStat("1111");
        H1F hi_dtheta_pos = new H1F("hi_dtheta_pos","hi_dtheta_pos", 100, -1, 1); 
        hi_dtheta_pos.setTitleX("#Delta #theta (deg)");
        hi_dtheta_pos.setTitleY("Counts");
        hi_dtheta_pos.setTitle("Positive Tracks");
        F1D f1_dtheta_pos = new F1D("f1_dtheta_pos","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dtheta_pos.setParameter(0, 0);
        f1_dtheta_pos.setParameter(1, 0);
        f1_dtheta_pos.setParameter(2, 1.0);
        f1_dtheta_pos.setLineWidth(2);
        f1_dtheta_pos.setLineColor(2);
        f1_dtheta_pos.setOptStat("1111");        
        H1F hi_dphi_pos = new H1F("hi_dphi_pos", "hi_dphi_pos", 100, -5.0, 5.0); 
        hi_dphi_pos.setTitleX("#Delta #phi (deg)");
        hi_dphi_pos.setTitleY("Counts");
        hi_dphi_pos.setTitle("Positive Tracks");
        F1D f1_dphi_pos = new F1D("f1_dphi_pos","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dphi_pos.setParameter(0, 0);
        f1_dphi_pos.setParameter(1, 0);
        f1_dphi_pos.setParameter(2, 1.0);
        f1_dphi_pos.setLineWidth(2);
        f1_dphi_pos.setLineColor(2);
        f1_dphi_pos.setOptStat("1111");        
        H1F hi_dvz_pos = new H1F("hi_dvz_pos", "hi_dvz_pos", 100, -10.0, 10.0);   
        hi_dvz_pos.setTitleX("#Delta Vz (cm)");
        hi_dvz_pos.setTitleY("Counts");
        hi_dvz_pos.setTitle("Positive Tracks");
        F1D f1_dvz_pos = new F1D("f1_dvz_pos","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dvz_pos.setParameter(0, 0);
        f1_dvz_pos.setParameter(1, 0);
        f1_dvz_pos.setParameter(2, 1.0);
        f1_dvz_pos.setLineWidth(2);
        f1_dvz_pos.setLineColor(2);
        f1_dvz_pos.setOptStat("1111");        
        H1F hi_dp_neg = new H1F("hi_dp_neg", "hi_dp_neg", 100, -0.1, 0.1); 
        hi_dp_neg.setTitleX("#Delta P/P");
        hi_dp_neg.setTitleY("Counts");
        hi_dp_neg.setTitle("Negative Tracks");
        F1D f1_dp_neg = new F1D("f1_dp_neg","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dp_neg.setParameter(0, 0);
        f1_dp_neg.setParameter(1, 0);
        f1_dp_neg.setParameter(2, 1.0);
        f1_dp_neg.setLineWidth(2);
        f1_dp_neg.setLineColor(2);
        f1_dp_neg.setOptStat("1111");
        H2F hi_dp_phi_neg = new H2F("hi_dp_phi_neg", "hi_dp_phi_neg", 100, -0.1, 0.1, 100, -180,180); 
        hi_dp_phi_neg.setTitleX("#Delta P/P");
        hi_dp_phi_neg.setTitleY("#phi (deg)");
        hi_dp_phi_neg.setTitle("Negative Tracks");
        H1F hi_dtheta_neg = new H1F("hi_dtheta_neg","hi_dtheta_neg", 100, -1, 1); 
        hi_dtheta_neg.setTitleX("#Delta #theta (deg)");
        hi_dtheta_neg.setTitleY("Counts");
        hi_dtheta_neg.setTitle("Negative Tracks");
        F1D f1_dtheta_neg = new F1D("f1_dtheta_neg","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dtheta_neg.setParameter(0, 0);
        f1_dtheta_neg.setParameter(1, 0);
        f1_dtheta_neg.setParameter(2, 1.0);
        f1_dtheta_neg.setLineWidth(2);
        f1_dtheta_neg.setLineColor(2);
        f1_dtheta_neg.setOptStat("1111");
        H1F hi_dphi_neg = new H1F("hi_dphi_neg", "hi_dphi_neg", 100, -5.0, 5.0); 
        hi_dphi_neg.setTitleX("#Delta #phi (deg)");
        hi_dphi_neg.setTitleY("Counts");
        hi_dphi_neg.setTitle("Negative Tracks");
        F1D f1_dphi_neg = new F1D("f1_dphi_neg","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dphi_neg.setParameter(0, 0);
        f1_dphi_neg.setParameter(1, 0);
        f1_dphi_neg.setParameter(2, 1.0);
        f1_dphi_neg.setLineWidth(2);
        f1_dphi_neg.setLineColor(2);
        f1_dphi_neg.setOptStat("1111");
        H1F hi_dvz_neg = new H1F("hi_dvz_neg", "hi_dvz_neg", 100, -10.0, 10.0);   
        hi_dvz_neg.setTitleX("#Delta Vz (cm)");
        hi_dvz_neg.setTitleY("Counts");
        hi_dvz_neg.setTitle("Negative Tracks");
        F1D f1_dvz_neg = new F1D("f1_dvz_neg","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dvz_neg.setParameter(0, 0);
        f1_dvz_neg.setParameter(1, 0);
        f1_dvz_neg.setParameter(2, 1.0);
        f1_dvz_neg.setLineWidth(2);
        f1_dvz_neg.setLineColor(2);
        f1_dvz_neg.setOptStat("1111");
        DataGroup mc = new DataGroup(4,2);
        mc.addDataSet(hi_dp_pos, 0);
        mc.addDataSet(f1_dp_pos, 0);
        mc.addDataSet(hi_dtheta_pos, 1);
        mc.addDataSet(f1_dtheta_pos, 1);
        mc.addDataSet(hi_dphi_pos, 2);
        mc.addDataSet(f1_dphi_pos, 2);
        mc.addDataSet(hi_dvz_pos, 3);
        mc.addDataSet(f1_dvz_pos, 3);
        mc.addDataSet(hi_dp_neg, 4);
        mc.addDataSet(f1_dp_neg, 4);
        mc.addDataSet(hi_dp_phi_neg, 0);
        mc.addDataSet(hi_dtheta_neg, 5);
        mc.addDataSet(f1_dtheta_neg, 5);
        mc.addDataSet(hi_dphi_neg, 6);
        mc.addDataSet(f1_dphi_neg, 6);
        mc.addDataSet(hi_dvz_neg, 7);
        mc.addDataSet(f1_dvz_neg, 7);
        this.getDataGroup().add(mc, 0);
        // Truth
        DataGroup dc_truth = new DataGroup(2, 3);
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            H1F hi_truth_res = new H1F("hi_truth_res_l" + layer, "Residual (cm) - Layer " + layer, "Counts", 200, -10, 10);
            hi_truth_res.setFillColor(4);
            dc_truth.addDataSet(hi_truth_res, layer - 1);
        }
        this.getDataGroup().add(dc_truth, 3);
        // vertex
        H2F hi_vxy_pos = new H2F("hi_vxy_pos","hi_vxy_pos",100,-15.,15.,100,-15.,15);
        hi_vxy_pos.setTitleX("Vx (cm)");
        hi_vxy_pos.setTitleY("Vy (cm)");
        H2F hi_vxy_neg = new H2F("hi_vxy_neg","hi_vxy_neg",100,-15.,15.,100,-15.,15);
        hi_vxy_neg.setTitleX("Vx (cm)");
        hi_vxy_neg.setTitleY("Vy (cm)"); 
        H2F hi_vz_vs_theta_pos = new H2F("hi_vz_vs_theta_pos","hi_vz_vs_theta_pos",100, 5.,40.,100,-15.,40);
        hi_vz_vs_theta_pos.setTitleX("#theta (deg)");
        hi_vz_vs_theta_pos.setTitleY("Vz (cm)");
        H2F hi_vz_vs_theta_neg = new H2F("hi_vz_vs_theta_neg","hi_vz_vs_theta_neg",100, 5.,40.,100,-15.,40);
        hi_vz_vs_theta_neg.setTitleX("#theta (deg)");
        hi_vz_vs_theta_neg.setTitleY("Vz (cm)");
        H2F hi_vz_vs_phi_pos = new H2F("hi_vz_vs_phi_pos","hi_vz_vs_phi_pos",200,-15.,15.,200,-180,180);
        hi_vz_vs_phi_pos.setTitleX("Vz (cm)");
        hi_vz_vs_phi_pos.setTitleY("#phi (deg)");
        H2F hi_vz_vs_phi_neg = new H2F("hi_vz_vs_phi_neg","hi_vz_vs_phi_neg",200,-15.,15.,200,-180,180);
        hi_vz_vs_phi_neg.setTitleX("Vz (cm)");
        hi_vz_vs_phi_neg.setTitleY("#phi (deg)");
        H1F hi_chi2_pos = new H1F("hi_chi2_pos", "hi_chi2_pos", 100, 0.0, 10.0);   
        hi_chi2_pos.setTitleX("#chi^2");
        hi_chi2_pos.setTitleY("Counts");
        hi_chi2_pos.setTitle("Negative Tracks");
        H1F hi_chi2_neg = new H1F("hi_chi2_neg", "hi_chi2_neg", 100, 0.0, 10.0);   
        hi_chi2_neg.setTitleX("#chi^2");
        hi_chi2_neg.setTitleY("Counts");
        hi_chi2_neg.setTitle("Negative Tracks");
        DataGroup vertex = new DataGroup(4,2);
        vertex.addDataSet(hi_vz_vs_theta_pos, 0);
        vertex.addDataSet(hi_vxy_pos, 1);
        vertex.addDataSet(hi_vz_vs_phi_pos, 2);
        vertex.addDataSet(hi_chi2_pos, 3);
        vertex.addDataSet(hi_vz_vs_theta_neg, 4);
        vertex.addDataSet(hi_vxy_neg, 5);
        vertex.addDataSet(hi_vz_vs_phi_neg, 6);
        vertex.addDataSet(hi_chi2_neg, 7);
        this.getDataGroup().add(vertex, 4);
        // Hits
        DataGroup dc_hits = new DataGroup(3, 4);
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            H1F hi_hit_res = new H1F("hi_hit_res_l" + layer, "Residual (cm) - Layer " + layer, "Counts", 200, -10, 10);
            hi_hit_res.setFillColor(4);
            H1F hi_hit_phi = new H1F("hi_hit_phi_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_hit_phi.setFillColor(4);
            H1F hi_hit_phi_cut = new H1F("hi_hit_phi_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_hit_phi_cut.setFillColor(3);
            dc_hits.addDataSet(hi_hit_res, layer - 1);
            dc_hits.addDataSet(hi_hit_phi, layer - 1 + FVT_Nlayers);
            dc_hits.addDataSet(hi_hit_phi_cut, layer - 1 + FVT_Nlayers);
        }
        this.getDataGroup().add(dc_hits, 1);
        // Clusters
        DataGroup dc_clusters = new DataGroup(3, 10);
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            H1F hi_cluster_res = new H1F("hi_cluster_res_l" + layer, "Residual (cm) - Layer " + layer, "Counts", 200, -10, 10);
            hi_cluster_res.setFillColor(4);
            F1D f1_cluster = new F1D("f1_cluster_l" + layer,"[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
            f1_cluster.setParameter(0, 0);
            f1_cluster.setParameter(1, 0);
            f1_cluster.setParameter(2, 1.0);
            f1_cluster.setLineWidth(2);
            f1_cluster.setLineColor(2);
            f1_cluster.setOptStat("1111");
            H1F hi_cluster_phi = new H1F("hi_cluster_phi_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_phi.setFillColor(4);
            H1F hi_cluster_phi_cut = new H1F("hi_cluster_phi_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_phi_cut.setFillColor(3);
            H1F hi_cluster_energy_cut = new H1F("hi_cluster_energy_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_energy_cut.setFillColor(2);
            H2F hi_cluster_phi_energy = new H2F("hi_cluster_phi_energy_l" + layer, 100, -180, 180, 100, 0, 500);
            hi_cluster_phi_energy.setTitleX("Strip angle - Layer " + layer);
            hi_cluster_phi_energy.setTitleY("Energy");
            H2F hi_cluster_phi_strip = new H2F("hi_cluster_phi_strip_l" + layer, 100, -180, 180, 100, 0, 1000);
            hi_cluster_phi_strip.setTitleX("Strip angle - Layer " + layer);
            hi_cluster_phi_strip.setTitleY("Strip - Layer " + layer);
            H2F hi_cluster_strip_strip = new H2F("hi_cluster_strip_strip_l" + layer, 100, 0, 1000, 100, 0, 1000);
            hi_cluster_strip_strip.setTitleX("Strip - Layer " + "3");
            hi_cluster_strip_strip.setTitleY("Strip - Layer " + layer);
            dc_clusters.addDataSet(hi_cluster_res, layer - 1);
            dc_clusters.addDataSet(f1_cluster, layer - 1);
            dc_clusters.addDataSet(hi_cluster_phi, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_cut, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_energy_cut, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_strip, layer - 1 + 2*FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_strip_strip, layer - 1 + 3*FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_energy, layer - 1 + 4*FVT_Nlayers);
        }
        this.getDataGroup().add(dc_clusters, 2);

    }

    @Override
    public void plotHistos() {

        this.getAnalysisCanvas().getCanvas("Monte Carlo").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Truth").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Truth").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Truth").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Vertex").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridY(false);
        this.getAnalysisCanvas().getCanvas("FMT hits").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("FMT hits").setGridX(false);
        this.getAnalysisCanvas().getCanvas("FMT hits").setGridY(false);
        this.getAnalysisCanvas().getCanvas("FMT phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("FMT phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("FMT phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Clusters").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Clusters").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Clusters").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster energy").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster energy").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster energy").setGridY(false);

        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(0);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dp_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dp_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(1);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dtheta_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dtheta_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dphi_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dphi_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(3);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dvz_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dvz_pos"),"same");
//        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(3);
//        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH2F("hi_dp_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(4);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dp_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dp_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(5);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dtheta_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dtheta_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(6);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dphi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dphi_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(7);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getH1F("hi_dvz_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(0).getF1D("f1_dvz_neg"),"same");
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.getAnalysisCanvas().getCanvas("Truth").cd(layer - 1);
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Truth").draw(this.getDataGroup().getItem(3).getH1F("hi_truth_res_l" + layer));
        }
        this.getAnalysisCanvas().getCanvas("Vertex").cd(0);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH1F("hi_chi2_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(1);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(2);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vxy_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(3);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_phi_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(4);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH1F("hi_chi2_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(5);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(6);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(6).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vxy_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(7);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_phi_neg"));
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.getAnalysisCanvas().getCanvas("FMT hits").cd(layer - 1);
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("FMT hits").draw(this.getDataGroup().getItem(1).getH1F("hi_hit_res_l" + layer));
            this.getAnalysisCanvas().getCanvas("FMT phi").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("FMT phi").draw(this.getDataGroup().getItem(1).getH1F("hi_hit_phi_l" + layer));
            this.getAnalysisCanvas().getCanvas("FMT phi").draw(this.getDataGroup().getItem(1).getH1F("hi_hit_phi_cut_l" + layer), "same");
        }
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.getAnalysisCanvas().getCanvas("Clusters").cd(layer - 1);
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Clusters").draw(this.getDataGroup().getItem(2).getH1F("hi_cluster_res_l" + layer));
            this.getAnalysisCanvas().getCanvas("Clusters").draw(this.getDataGroup().getItem(2).getF1D("f1_cluster_l" + layer),"same");
            this.getAnalysisCanvas().getCanvas("Cluster phi").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(2).getH1F("hi_cluster_energy_cut_l" + layer), "same");
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_cut_l" + layer), "same");
            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").cd(layer - 1);
//            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").getPad(layer - 1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").draw(this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_strip_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").getPad(layer - 1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").draw(this.getDataGroup().getItem(2).getH2F("hi_cluster_strip_strip_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster energy").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster energy").draw(this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_energy_l" + layer));
        }
    }

    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun = null;
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        DataBank recTrajEB = null;
        DataBank fmtHits = null;
        DataBank fmtClusters = null;
        DataBank fmtADC = null;
        DataBank mcTrue = null;
        if (event.hasBank("RUN::config")) {
            recRun = event.getBank("RUN::config");
        }
        if (event.hasBank("REC::Particle")) {
            recBankEB = event.getBank("REC::Particle");
        }
        if (event.hasBank("REC::Event")) {
            recEvenEB = event.getBank("REC::Event");
        }
        if (event.hasBank("REC::Traj")) {
            recTrajEB = event.getBank("REC::Traj");
        }
        if (event.hasBank("FMT::Hits")) {
            fmtHits = event.getBank("FMT::Hits");
        }
        if (event.hasBank("FMT::Clusters")) {
            fmtClusters = event.getBank("FMT::Clusters");
        }
        if (event.hasBank("FMT::adc")) {
            fmtADC = event.getBank("FMT::adc");
        }
        if (event.hasBank("MC::True")) {
            mcTrue = event.getBank("MC::True");
        }
        if(recRun == null) return;
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        if(run==0) return;
        IndexedTable rfConfig = this.getCcdb().getConstants(run, "/calibration/eb/rf/config");
        if(this.rfPeriod!=rfConfig.getDoubleValue("clock", 1,1,1)) {
            this.rfPeriod = rfConfig.getDoubleValue("clock", 1,1,1);
            this.resetEventListener();
        }
        // Decoding Trigger Bits
        boolean[] trigger_bits = new boolean[32];
        if (event.hasBank("RUN::config")) {
            DataBank bank = event.getBank("RUN::config");
            long TriggerWord = bank.getLong("trigger", 0) & 0xFFFFFFFF;
            for (int i = 31; i >= 0; i--) {
                trigger_bits[i] = (TriggerWord & (1 << i)) != 0;
            }
        }
        // get event start time
        double startTime = -1000;
        double rfTime = -1000;
        if (recEvenEB != null) {
            startTime = recEvenEB.getFloat("startTime", 0);
            rfTime = recEvenEB.getFloat("RFTime", 0);
        }
        // get trigger particle
        int trigger = 0;
        if (recBankEB != null) {
            trigger = recBankEB.getInt("pid", 0);
        }
        
        Particle partGenNeg = null;
        Particle partGenPos = null;
        Particle partRecNeg = null;
        Particle partRecPos = null;        
        if(event.hasBank("FMT::Tracks")==true){
            DataBank  bank = event.getBank("FMT::Tracks");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
                int pidCode=0;
                if(bank.getByte("q", loop)==-1) pidCode = 11;
                else                            pidCode = 211;
                Particle recParticle = new Particle(
                                          pidCode,
                                          bank.getFloat("p0_x", loop),
                                          bank.getFloat("p0_y", loop),
                                          bank.getFloat("p0_z", loop),
                                          bank.getFloat("Vtx0_x", loop),
                                          bank.getFloat("Vtx0_y", loop),
                                          bank.getFloat("Vtx0_z", loop));
                if(partRecNeg==null && recParticle.charge()<0) {
                    partRecNeg=recParticle;
                    this.getDataGroup().getItem(4).getH2F("hi_vxy_neg").fill(recParticle.vx(),recParticle.vy());
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_neg").fill(recParticle.vz(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_phi_neg").fill(recParticle.vz(),Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(4).getH1F("hi_chi2_neg").fill(bank.getFloat("chi2", loop));
                }
                if(partRecPos==null && recParticle.charge()>0) {
                    partRecPos=recParticle;
                    this.getDataGroup().getItem(4).getH2F("hi_vxy_pos").fill(recParticle.vx(),recParticle.vy());
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_pos").fill(recParticle.vz(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_phi_pos").fill(recParticle.vz(),Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(4).getH1F("hi_chi2_pos").fill(bank.getFloat("chi2", loop));
                }
            }
        }
        if(event.hasBank("MC::Particle")==true){
            DataBank genBank = event.getBank("MC::Particle");
            int nrows = genBank.rows();
            for(int loop = 0; loop < nrows; loop++) {   
                Particle genPart = new Particle(
                                              genBank.getInt("pid", loop),
                                              genBank.getFloat("px", loop),
                                              genBank.getFloat("py", loop),
                                              genBank.getFloat("pz", loop),
                                              genBank.getFloat("vx", loop),
                                              genBank.getFloat("vy", loop),
                                              genBank.getFloat("vz", loop));
                if(genPart.charge()==-1  && partGenNeg==null && genPart.theta()!=0) partGenNeg=genPart;
                if(genPart.charge()==1   && partGenPos==null) partGenPos=genPart;
            }
            if(partGenNeg != null && partRecNeg != null) {
                this.getDataGroup().getItem(0).getH1F("hi_dp_neg").fill((partRecNeg.p()-partGenNeg.p())/partGenNeg.p());
                this.getDataGroup().getItem(0).getH2F("hi_dp_phi_neg").fill((partRecNeg.p()-partGenNeg.p())/partGenNeg.p(),Math.toDegrees(partRecNeg.phi()));
                this.getDataGroup().getItem(0).getH1F("hi_dtheta_neg").fill(Math.toDegrees(partRecNeg.theta()-partGenNeg.theta()));
                this.getDataGroup().getItem(0).getH1F("hi_dphi_neg").fill(Math.toDegrees(partRecNeg.phi()-partGenNeg.phi()));
                this.getDataGroup().getItem(0).getH1F("hi_dvz_neg").fill(partRecNeg.vz()-partGenNeg.vz());
            }
            if(partGenPos != null && partRecPos != null) {
                this.getDataGroup().getItem(0).getH1F("hi_dp_pos").fill((partRecPos.p()-partGenPos.p())/partGenPos.p());
                this.getDataGroup().getItem(0).getH1F("hi_dtheta_pos").fill(Math.toDegrees(partRecPos.theta()-partGenPos.theta()));
                this.getDataGroup().getItem(0).getH1F("hi_dphi_pos").fill(Math.toDegrees(partRecPos.phi()-partGenPos.phi()));
                this.getDataGroup().getItem(0).getH1F("hi_dvz_pos").fill(partRecPos.vz()-partGenPos.vz());
            }
            
        }
        if(recTrajEB != null) {
            for (int loop = 0; loop < recTrajEB.rows(); loop++) {
                int detId = recTrajEB.getInt("detector", loop);
                int layer = recTrajEB.getByte("layer", loop);
                if (detId==DetectorType.FMT.getDetectorId() && layer==1) {
                    double x = recTrajEB.getFloat("x", loop);
                    double y = recTrajEB.getFloat("y", loop);
                    double z = recTrajEB.getFloat("z", loop);
                    double angle = Math.toDegrees(Math.acos(z/Math.sqrt(x*x+y*y+z*z)));
                    if(angle>7 && angle<25) {
                        ntraj++;
                        if(fmtHits!=null) nfmt++;
                    }
                }
            }
        }

        if (fmtHits != null && fmtClusters != null && recTrajEB != null) {
//            System.out.println("Updating AAA");
            int[] matchedLayers = {0, 0, 0, 0, 0, 0};
            List<Integer> matchedStrips = new ArrayList<Integer>();
            List<Integer> matchedClusters = new ArrayList<Integer>();
            for (int loop = 0; loop < recTrajEB.rows(); loop++) {
                int detId = recTrajEB.getByte("detector", loop);
                int layer = recTrajEB.getByte("layer", loop);
                if (detId==DetectorType.FMT.getDetectorId() && layer>= 1 && layer <= FVT_Nlayers) {
                    double x0 = recTrajEB.getFloat("x", loop);
                    double y0 = recTrajEB.getFloat("y", loop);
                    double x = x0;
                    double y = y0;
                    double z = recTrajEB.getFloat("z", loop);
                    double phiRef = 0;
                    if(mcTrue!=null && fmtADC!=null) {
                        for(int i=0; i<mcTrue.rows(); i++) {
                            if(layer==fmtADC.getByte("layer", i)) {
                                this.getDataGroup().getItem(3).getH1F("hi_truth_res_l" + layer).fill(x0 - mcTrue.getFloat("avgX", i)/10);
                                this.getDataGroup().getItem(3).getH1F("hi_truth_res_l" + (layer+3)).fill(y0 - mcTrue.getFloat("avgY", i)/10);
                            }
                        }
                    }
                    if (event.hasBank("MC::Particle")) {
                        phiRef = FVT_Alpha[layer - 1];
                    } else {
                        phiRef = MY_Alpha[layer - 1];
                    }
                    for (int i = 0; i < fmtHits.rows(); i++) {
                        if (layer == fmtHits.getByte("layer", i)) {
                            int strip = fmtHits.getInt("strip", i);
                            double xLoc = x * Math.cos(Math.toRadians(phiRef)) + y * Math.sin(Math.toRadians(phiRef));
                            double yLoc = y * Math.cos(Math.toRadians(phiRef)) - x * Math.sin(Math.toRadians(phiRef));
                            double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
//                            if (Math.abs(yLoc - FVT_stripsYlocref[strip]) < 5) {
//                                matchedLayers[detId - 1]++;
//                            }
                            if(layer==3) matchedStrips.add(strip);
                            this.getDataGroup().getItem(1).getH1F("hi_hit_res_l" + layer).fill(yLoc - FVT_stripsYlocref[strip]);
                            if (phi[0] > -9999 && strip<=512) {
                                this.getDataGroup().getItem(1).getH1F("hi_hit_phi_l" + layer).fill(phi[0]);
                            }
                            if (phi[1] > -9999 && strip>512) {
                                this.getDataGroup().getItem(1).getH1F("hi_hit_phi_l" + layer).fill(phi[1]);
                            }
                        }
                    }
                    for (int i = 0; i < fmtClusters.rows(); i++) {
                        if (layer == fmtClusters.getByte("layer", i)) {
                            int strip = fmtClusters.getInt("seedStrip", i);
                            double energy = fmtClusters.getFloat("energy", i);
                            double xLoc = x * Math.cos(Math.toRadians(phiRef)) + y * Math.sin(Math.toRadians(phiRef));
                            double yLoc = y * Math.cos(Math.toRadians(phiRef)) - x * Math.sin(Math.toRadians(phiRef));
                            double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
//                            if (Math.abs(yLoc - FVT_stripsYlocref[strip]) < 5) {
//                                matchedLayers[detId - 1]++;
//                            }
                            if(layer==3) matchedClusters.add(strip);
//                            this.getDataGroup().getItem(2).getH1F("hi_cluster_res_l" + (layer+3)).fill(yLoc - FVT_stripsYlocref[strip]);
                            this.getDataGroup().getItem(2).getH1F("hi_cluster_res_l" + (layer+3)).fill(yLoc - fmtClusters.getFloat("centroid", i));
                            this.getDataGroup().getItem(2).getH1F("hi_cluster_res_l" + layer).fill(fmtClusters.getFloat("residual", i));
                            if (phi[0] > -9999 && strip<=512) {
                                this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_l" + layer).fill(phi[0]);
                                this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_energy_l" + layer).fill(phi[0],energy);
                                if(energy>100) this.getDataGroup().getItem(2).getH1F("hi_cluster_energy_cut_l" + layer).fill(phi[0]);
                                this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_strip_l" + layer).fill(phi[0],strip);
                            }
                            if (phi[1] > -9999 && strip>512) {
                                this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_l" + layer).fill(phi[1]);
                                this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_energy_l" + layer).fill(phi[1],energy);
                                if(energy>100) this.getDataGroup().getItem(2).getH1F("hi_cluster_energy_cut_l" + layer).fill(phi[1]);
                                this.getDataGroup().getItem(2).getH2F("hi_cluster_phi_strip_l" + layer).fill(phi[1],strip);
                            }
                        }
                    }
                }
            }
            for (int loop = 0; loop < recTrajEB.rows(); loop++) {
                int detId = recTrajEB.getByte("detector", loop);
                int layer = recTrajEB.getByte("layer", loop);
                if (detId==DetectorType.FMT.getDetectorId() && layer>= 1 && layer <= FVT_Nlayers) {
                    double x = recTrajEB.getFloat("x", loop);
                    double y = recTrajEB.getFloat("y", loop);
                    double z = recTrajEB.getFloat("z", loop);
//                    if (matchedLayers[0] > 0 && matchedLayers[1] > 0 && matchedLayers[3] > 0 && matchedLayers[4] > 0) {
                    for (int i = 0; i < fmtHits.rows(); i++) {
                        if (layer == fmtHits.getByte("layer", i)) {
                            int strip = fmtHits.getInt("strip", i);
                            boolean match = false;
                            for(int strip3 : matchedStrips) if(layer==6 && Math.abs(strip-strip3)<30) match=true;
                            if(match) {
                                double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
                                if (phi[0] > -9999 && strip<=512) {
                                    this.getDataGroup().getItem(1).getH1F("hi_hit_phi_cut_l" + layer).fill(phi[0]);
                                }
                                if (phi[1] > -9999 && strip>512) {
                                    this.getDataGroup().getItem(1).getH1F("hi_hit_phi_cut_l" + layer).fill(phi[1]);
                                }
                            }
                        }
                    }
                    for (int i = 0; i < fmtClusters.rows(); i++) {
                        if (layer == fmtClusters.getByte("layer", i)) {
                            int strip = fmtClusters.getInt("seedStrip", i);
                            double energy = fmtClusters.getFloat("energy", i);
                            boolean match = false;
                            for(int strip3 : matchedClusters) {
                                if(layer==6 && Math.abs(strip-strip3)<30) match=true;
                                this.getDataGroup().getItem(2).getH2F("hi_cluster_strip_strip_l" + layer).fill(strip3,strip);
                            }
                            if(match) {
                                double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
                                if (phi[0] > -9999 && strip<=512 && energy>100) {
                                    this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_cut_l" + layer).fill(phi[0]);
                                }
                                if (phi[1] > -9999 && strip>512 && energy>100) {
                                    this.getDataGroup().getItem(2).getH1F("hi_cluster_phi_cut_l" + layer).fill(phi[1]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating FTOF");
        this.analyze();
        //fitting negative tracks vertex
//        ParallelSliceFitter fitter = new ParallelSliceFitter(this.getDataGroup().getItem(4).getH2F("hi_rf_paddle_2"));
//        fitter.fitSlicesX();
//        GraphErrors meanY = fitter.getMeanSlices();
//        this.getDataGroup().getItem(4).getGraph("hi_rf_offsets_2").copy(meanY);       
    }

    @Override
    public void analyze() {
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dp_pos"),     this.getDataGroup().getItem(0).getF1D("f1_dp_pos"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dtheta_pos"), this.getDataGroup().getItem(0).getF1D("f1_dtheta_pos"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dphi_pos"),   this.getDataGroup().getItem(0).getF1D("f1_dphi_pos"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dvz_pos"),    this.getDataGroup().getItem(0).getF1D("f1_dvz_pos"));     
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dp_neg"),     this.getDataGroup().getItem(0).getF1D("f1_dp_neg"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dtheta_neg"), this.getDataGroup().getItem(0).getF1D("f1_dtheta_neg"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dphi_neg"),   this.getDataGroup().getItem(0).getF1D("f1_dphi_neg"));
        this.fitMC(this.getDataGroup().getItem(0).getH1F("hi_dvz_neg"),    this.getDataGroup().getItem(0).getF1D("f1_dvz_neg"));     
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.fitMC(this.getDataGroup().getItem(2).getH1F("hi_cluster_res_l" + layer), this.getDataGroup().getItem(2).getF1D("f1_cluster_l" + layer));
        }     
   }

    public void fitMC(H1F himc, F1D f1mc) {
        double mean  = himc.getDataX(himc.getMaximumBin());
        double amp   = himc.getBinContent(himc.getMaximumBin());
        double sigma = himc.getRMS();
        f1mc.setParameter(0, amp);
        f1mc.setParameter(1, mean);
        f1mc.setParameter(2, sigma);
        f1mc.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1mc, himc, "Q"); //No options uses error for sigma 
        sigma = Math.abs(f1mc.getParameter(2));  
        f1mc.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1mc, himc, "Q"); //No options uses error for sigma 
    }

    public void LoadGeometry() {

        
        for(int i=0; i<FVT_Nlayers; i++) FVT_Zlayer[i]+=2.5;
        
        int FVT_Sidestrips = (FVT_Nstrips - 2 * FVT_Halfstrips) / 2;
        double FVT_YCentral = (double) FVT_Halfstrips * FVT_Pitch / 2.;
        double FVT_Rmax = FVT_Pitch * (FVT_Halfstrips + 2 * FVT_Sidestrips) / 2.;
        double FVT_SigmaS = FVT_Pitch / Math.sqrt(12);
        double[][] FVT_stripsXloc = new double[FVT_Nstrips][2];
        double[][] FVT_stripsYloc = new double[FVT_Nstrips][2];
        FVT_stripsXlocref = new double[FVT_Nstrips];
        FVT_stripsYlocref = new double[FVT_Nstrips];
        FVT_stripsX = new double[FVT_Nlayers][FVT_Nstrips][2];
        FVT_stripsY = new double[FVT_Nlayers][FVT_Nstrips][2];
        FVT_stripslength = new double[FVT_Nstrips];

        for (int i = 0; i < FVT_Nstrips; i++) {
            //Give the Y of the middle of the strip
            if (i<512){
                    FVT_stripsYloc[i][0]=-FVT_Rmax+(511-i+0.5)*FVT_Pitch;
                    FVT_stripsYloc[i][1]=-FVT_Rmax+(511-i+0.5)*FVT_Pitch;
            } else {
                    FVT_stripsYloc[i][0]=FVT_Rmax-(1023-i+0.5)*FVT_Pitch;
                    FVT_stripsYloc[i][1]=FVT_Rmax-(1023-i+0.5)*FVT_Pitch;
            }
            FVT_stripsYlocref[i] = FVT_stripsYloc[i][0];

            int localRegion = getLocalRegion(i);
            switch (localRegion) {
                case 2:
                case 4:
                    FVT_stripslength[i] = 2 * FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][0] = -FVT_stripslength[i] / 2.;
                    FVT_stripsXloc[i][1] = FVT_stripslength[i] / 2.;
                    FVT_stripsXlocref[i] = 0;
                    break;
                case 1:
                    FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][1] = 0;
                    FVT_stripsXloc[i][0] = -FVT_stripslength[i];
                    FVT_stripsXlocref[i] = -FVT_stripslength[i] / 2;
                    if (Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole < 1) {
                        FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax)) - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][1] = -FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][0] = -FVT_stripslength[i];
                        FVT_stripsXlocref[i] = -FVT_stripslength[i] / 2 - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                    }
                    break;
                case 3:
                    FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][0] = 0;
                    FVT_stripsXloc[i][1] = FVT_stripslength[i];
                    FVT_stripsXlocref[i] = FVT_stripslength[i] / 2;
                    if (Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole < 1) {
                        FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax)) - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][0] = FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][1] = FVT_stripslength[i];
                        FVT_stripsXlocref[i] = FVT_stripslength[i] / 2 + FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                    }
                    break;
            }
            for (int j = 0; j < FVT_Nlayers; j++) { //X sign flip
                FVT_stripsX[j][i][0] = -(FVT_stripsXloc[i][0] * Math.cos(FVT_Alpha[j]) + FVT_stripsYloc[i][0] * Math.sin(FVT_Alpha[j]));
                FVT_stripsY[j][i][0] = -FVT_stripsXloc[i][0] * Math.sin(FVT_Alpha[j]) + FVT_stripsYloc[i][0] * Math.cos(FVT_Alpha[j]);
                FVT_stripsX[j][i][1] = -(FVT_stripsXloc[i][1] * Math.cos(FVT_Alpha[j]) + FVT_stripsYloc[i][1] * Math.sin(FVT_Alpha[j]));
                FVT_stripsY[j][i][1] = -FVT_stripsXloc[i][1] * Math.sin(FVT_Alpha[j]) + FVT_stripsYloc[i][1] * Math.cos(FVT_Alpha[j]);
            }

            //System.out.println(Constants.getLocalRegion(i)+" strip-1 = "+i+" x' "+FVT_stripsXloc[i][1]+" y' "+FVT_stripsYloc[i][1]+" length "+FVT_stripslength[i]+" FVT_Beamhole "+FVT_Beamhole);
        }

    }

    private int getLocalRegion(int i) {
		// To represent the geometry we divide the barrel micromega disk into 3 regions according to the strip numbering system.
        // Here i = strip_number -1;
        // Region 1 is the region in the negative x part of inner region: the strips range is from   1 to 320  (   0 <= i < 320)
        // Region 2 is the region in the negative y part of outer region: the strips range is from 321 to 512  ( 320 <= i < 512)
        // Region 3 is the region in the positive x part of inner region: the strips range is from 513 to 832  ( 512 <= i < 832)
        // Region 4 is the region in the positive y part of outer region: the strips range is from 833 to 1024 ( 832 <= i < 1024)

        int region = 0;
        if (i >= 0 && i < 320) {
            region = 1;
        }
        if (i >= 320 && i < 512) {
            region = 2;
        }
        if (i >= 512 && i < 832) {
            region = 3;
        }
        if (i >= 832 && i < 1024) {
            region = 4;
        }

        return region;
    }

    private double[] getPhi(double x, double y, double yloc) {
        double t1 = -9999;
        double t2 = -9999;
        if (y == -yloc) {
            t1 = y / x;
            t2 = y / x;
        } else if (x * x > (yloc * yloc - y * y)) {
            t1 = (-x + Math.sqrt(x * x - (yloc * yloc - y * y))) / (yloc + y);
            t2 = (-x - Math.sqrt(x * x - (yloc * yloc - y * y))) / (yloc + y);
//            if(Math.abs(yloc-y*(1-t1*t1)/(1+t1*t1))-x*2*t1/(1+t1*t1)<0.1 || true) t=t1;
//            else t=t2;
        }
        double phi[] = new double[2];
        phi[0]=2 * Math.toDegrees(Math.atan(t1));
        phi[1]=2 * Math.toDegrees(Math.atan(t2));
        if(t1==-9999) phi[0]=-9999;
        if(t2==-9999) phi[1]=-9999;
        return phi;
    }
}
