/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.detector.base.DetectorType;

/**
 *
 * @author ziegler
 */
public class CVTmonitor extends AnalysisMonitor {
    

    public CVTmonitor(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames("Monte Carlo","Monte Carlo 2D","Vertex","Positive Tracks","Negative Tracks","Trajectories");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getAnalysisCanvas().getCanvas("Trajectories").divide(4,2);
        this.getAnalysisCanvas().getCanvas("Trajectories").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Trajectories").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").divide(4,2);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").divide(4,2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Vertex").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridY(false);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY(" hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // negative tracks
        H1F hi_p_neg = new H1F("hi_p_neg", "hi_p_neg", 100, 0.0, 3.0);     
        hi_p_neg.setTitleX("p (GeV)");
        hi_p_neg.setTitleY("Counts");
        H1F hi_theta_neg = new H1F("hi_theta_neg", "hi_theta_neg", 100, 30.0, 150.0); 
        hi_theta_neg.setTitleX("#theta (deg)");
        hi_theta_neg.setTitleY("Counts");
        H1F hi_phi_neg = new H1F("hi_phi_neg", "hi_phi_neg", 100, -180.0, 180.0);   
        hi_phi_neg.setTitleX("#phi (deg)");
        hi_phi_neg.setTitleY("Counts");
        H1F hi_chi2_neg = new H1F("hi_chi2_neg", "hi_chi2_neg", 100, 0, 50);   
        hi_chi2_neg.setTitleX("#chi2");
        hi_chi2_neg.setTitleY("Counts");
        H1F hi_vz_neg = new H1F("hi_vz_neg", "hi_vz_neg", 100, -10.0, 10.0);   
        hi_vz_neg.setTitleX("Vz (cm)");
        hi_vz_neg.setTitleY("Counts");
        H2F hi_theta_p_neg = new H2F("hi_theta_p_neg", "hi_theta_p_neg", 100, 0.0, 3.0, 100, 30.0, 150.0); 
        hi_theta_p_neg.setTitleX("p (GeV)");
        hi_theta_p_neg.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_neg = new H2F("hi_theta_phi_neg", "hi_theta_phi_neg", 100, -180.0, 180.0, 100, 30.0, 150.0); 
        hi_theta_phi_neg.setTitleX("#phi (deg)");
        hi_theta_phi_neg.setTitleY("#theta (deg)");        
        H1F hi_ndf_neg = new H1F("hi_ndf_neg", "hi_ndf_neg", 15, 0, 15);   
        hi_ndf_neg.setTitleX("NDF");
        hi_ndf_neg.setTitleY("Counts");
        DataGroup dg_neg = new DataGroup(3,2);
        dg_neg.addDataSet(hi_p_neg,         0);
        dg_neg.addDataSet(hi_theta_neg,     1);
        dg_neg.addDataSet(hi_phi_neg,       2);
        dg_neg.addDataSet(hi_chi2_neg,      3);
        dg_neg.addDataSet(hi_vz_neg,        4);
        dg_neg.addDataSet(hi_theta_p_neg,   5);
        dg_neg.addDataSet(hi_theta_phi_neg, 6);
        dg_neg.addDataSet(hi_ndf_neg,       7);
        this.getDataGroup().add(dg_neg, 1);
        // positive trakcs
        H1F hi_p_pos = new H1F("hi_p_pos", "hi_p_pos", 100, 0.0, 3.0);     
        hi_p_pos.setTitleX("p (GeV)");
        hi_p_pos.setTitleY("Counts");
        H1F hi_theta_pos = new H1F("hi_theta_pos", "hi_theta_pos", 100, 30.0, 150.0); 
        hi_theta_pos.setTitleX("#theta (deg)");
        hi_theta_pos.setTitleY("Counts");
        H1F hi_phi_pos = new H1F("hi_phi_pos", "hi_phi_pos", 100, -180.0, 180.0);   
        hi_phi_pos.setTitleX("#phi (deg)");
        hi_phi_pos.setTitleY("Counts");
        H1F hi_chi2_pos = new H1F("hi_chi2_pos", "hi_chi2_pos", 100, 0, 50);   
        hi_chi2_pos.setTitleX("#chi2");
        hi_chi2_pos.setTitleY("Counts");
        H1F hi_vz_pos = new H1F("hi_vz_pos", "hi_vz_pos", 100, -10.0, 10.0);   
        hi_vz_pos.setTitleX("Vz (cm)");
        hi_vz_pos.setTitleY("Counts");
        H2F hi_theta_p_pos = new H2F("hi_theta_p_pos", "hi_theta_p_pos", 100, 0.0, 3.0, 100, 30.0, 150.0); 
        hi_theta_p_pos.setTitleX("p (GeV)");
        hi_theta_p_pos.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_pos = new H2F("hi_theta_phi_pos", "hi_theta_phi_pos", 100, -180.0, 180.0, 100, 30.0, 150.0); 
        hi_theta_phi_pos.setTitleX("#phi (deg)");
        hi_theta_phi_pos.setTitleY("#theta (deg)");  
        H1F hi_ndf_pos = new H1F("hi_ndf_pos", "hi_ndf_pos", 15, 0, 15);   
        hi_ndf_pos.setTitleX("NDF");
        hi_ndf_pos.setTitleY("Counts");
        // mc positive trakcs
        H1F hi_mc_p_pos = new H1F("hi_mc_p_pos", "hi_mc_p_pos", 100, 0.0, 3.0);     
        hi_mc_p_pos.setTitleX("p (GeV)");
        hi_mc_p_pos.setTitleY("Counts");
        hi_mc_p_pos.setLineColor(2);
        H1F hi_mc_theta_pos = new H1F("hi_mc_theta_pos", "hi_mc_theta_pos", 100, 30.0, 150.0); 
        hi_mc_theta_pos.setTitleX("#theta (deg)");
        hi_mc_theta_pos.setTitleY("Counts");
        hi_mc_theta_pos.setLineColor(2);
        H1F hi_mc_phi_pos = new H1F("hi_mc_phi_pos", "hi_mc_phi_pos", 100, -180.0, 180.0);   
        hi_mc_phi_pos.setTitleX("#phi (deg)");
        hi_mc_phi_pos.setTitleY("Counts");
        hi_mc_phi_pos.setLineColor(2);
        H1F hi_mc_vz_pos = new H1F("hi_mc_vz_pos", "hi_mc_vz_pos", 100, -10.0, 10.0);   
        hi_mc_vz_pos.setTitleX("Vz (cm)");
        hi_mc_vz_pos.setTitleY("Counts");
        hi_mc_vz_pos.setLineColor(2);
        DataGroup dg_pos = new DataGroup(3,2);
        dg_pos.addDataSet(hi_p_pos,         0);
        dg_pos.addDataSet(hi_mc_p_pos,      0);
        dg_pos.addDataSet(hi_theta_pos,     1);
        dg_pos.addDataSet(hi_mc_theta_pos,  1);
        dg_pos.addDataSet(hi_phi_pos,       2);
        dg_pos.addDataSet(hi_mc_phi_pos,    2);
        dg_pos.addDataSet(hi_chi2_pos,      3);
        dg_pos.addDataSet(hi_vz_pos,        4);
        dg_pos.addDataSet(hi_mc_vz_pos,     4);
        dg_pos.addDataSet(hi_theta_p_pos,   5);
        dg_pos.addDataSet(hi_theta_phi_pos, 6);
        dg_pos.addDataSet(hi_ndf_pos,       7);
        this.getDataGroup().add(dg_pos, 2);
        // mc comparison
        H1F hi_dp_pos = new H1F("hi_dp_pos", "hi_dp_pos", 100, -0.4, 0.4); 
        hi_dp_pos.setTitleX("#Delta P / P ");
        hi_dp_pos.setTitleY("Counts");
        hi_dp_pos.setTitle("Positive Tracks");
        F1D f1_dp_pos = new F1D("f1_dp_pos","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dp_pos.setParameter(0, 0);
        f1_dp_pos.setParameter(1, 0);
        f1_dp_pos.setParameter(2, 1.0);
        f1_dp_pos.setLineWidth(2);
        f1_dp_pos.setLineColor(2);
        f1_dp_pos.setOptStat("1111");
        H1F hi_dtheta_pos = new H1F("hi_dtheta_pos","hi_dtheta_pos", 100, -4.0, 4.0); 
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
        H1F hi_dphi_pos = new H1F("hi_dphi_pos", "hi_dphi_pos", 100, -4.0, 4.0); 
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
        H1F hi_dvz_pos = new H1F("hi_dvz_pos", "hi_dvz_pos", 100, -5.0, 5.0);   
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
        H1F hi_dp_neg = new H1F("hi_dp_neg", "hi_dp_neg", 100, -0.4, 0.4); 
        hi_dp_neg.setTitleX("#Delta P / P ");
        hi_dp_neg.setTitleY("Counts");
        hi_dp_neg.setTitle("Negative Tracks");
        F1D f1_dp_neg = new F1D("f1_dp_neg","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        f1_dp_neg.setParameter(0, 0);
        f1_dp_neg.setParameter(1, 0);
        f1_dp_neg.setParameter(2, 1.0);
        f1_dp_neg.setLineWidth(2);
        f1_dp_neg.setLineColor(2);
        f1_dp_neg.setOptStat("1111");
        H1F hi_dtheta_neg = new H1F("hi_dtheta_neg","hi_dtheta_neg", 100, -4.0, 4.0); 
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
        H1F hi_dphi_neg = new H1F("hi_dphi_neg", "hi_dphi_neg", 100, -4.0, 4.0); 
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
        H1F hi_dvz_neg = new H1F("hi_dvz_neg", "hi_dvz_neg", 100, -5.0, 5.0);   
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
        mc.addDataSet(hi_dtheta_neg, 5);
        mc.addDataSet(f1_dtheta_neg, 5);
        mc.addDataSet(hi_dphi_neg, 6);
        mc.addDataSet(f1_dphi_neg, 6);
        mc.addDataSet(hi_dvz_neg, 7);
        mc.addDataSet(f1_dvz_neg, 7);
        this.getDataGroup().add(mc, 3);
        // MC-2D
        H2F hi_dp_p_pos = new H2F("hi_dp_p_pos","hi_dp_p_pos",100, 0., 3., 100, -0.4, 0.4);
        hi_dp_p_pos.setTitleX("p (GeV)");
        hi_dp_p_pos.setTitleY("#DeltaP/P ");
        H2F hi_dp_theta_pos = new H2F("hi_dp_theta_pos","hi_dp_theta_pos",100, 35, 155, 100, -0.4, 0.4);
        hi_dp_theta_pos.setTitleX("#theta (deg)");
        hi_dp_theta_pos.setTitleY("#DeltaP/P ");
        H2F hi_dp_phi_pos = new H2F("hi_dp_phi_pos","hi_dp_phi_pos",100, -180, 180, 100, -0.4, 0.4);
        hi_dp_phi_pos.setTitleX("#phi (deg)");
        hi_dp_phi_pos.setTitleY("#DeltaP/P ");
        H2F hi_dp_z_pos = new H2F("hi_dp_z_pos","hi_dp_z_pos",100, -10, 10, 100, -0.4, 0.4);
        hi_dp_z_pos.setTitleX("z (cm)");
        hi_dp_z_pos.setTitleY("#DeltaP/P ");
        H2F hi_dp_p_neg = new H2F("hi_dp_p_neg","hi_dp_p_neg",100, 0., 3., 100, -0.4, 0.4);
        hi_dp_p_neg.setTitleX("p (GeV)");
        hi_dp_p_neg.setTitleY("#DeltaP/P ");
        H2F hi_dp_theta_neg = new H2F("hi_dp_theta_neg","hi_dp_theta_neg",100, 35, 155, 100, -0.4, 0.4);
        hi_dp_theta_neg.setTitleX("#theta (deg)");
        hi_dp_theta_neg.setTitleY("#DeltaP/P ");
        H2F hi_dp_phi_neg = new H2F("hi_dp_phi_neg","hi_dp_phi_neg",100, -180, 180, 100, -0.4, 0.4);
        hi_dp_phi_neg.setTitleX("#phi (deg)");
        hi_dp_phi_neg.setTitleY("#DeltaP/P ");
        H2F hi_dp_z_neg = new H2F("hi_dp_z_neg","hi_dp_z_neg",100, -10, 10, 100, -0.4, 0.4);
        hi_dp_z_neg.setTitleX("z (cm)");
        hi_dp_z_neg.setTitleY("#DeltaP/P ");
        DataGroup mc2d = new DataGroup(4,2);
        mc2d.addDataSet(hi_dp_p_pos,     0);
        mc2d.addDataSet(hi_dp_theta_pos, 1);
        mc2d.addDataSet(hi_dp_phi_pos,   2);
        mc2d.addDataSet(hi_dp_z_pos,     3);
        mc2d.addDataSet(hi_dp_p_neg,     4);
        mc2d.addDataSet(hi_dp_theta_neg, 5);
        mc2d.addDataSet(hi_dp_phi_neg,   6);
        mc2d.addDataSet(hi_dp_z_neg,     7);
        this.getDataGroup().add(mc2d, 4);
        // vertex
        H2F hi_vxy_pos = new H2F("hi_vxy_pos","hi_vxy_pos",100,-1.,1.,100,-1.,1);
        hi_vxy_pos.setTitleX("Vx (cm)");
        hi_vxy_pos.setTitleY("Vy (cm)");
        H2F hi_vxy_neg = new H2F("hi_vxy_neg","hi_vxy_neg",100,-0.2,0.2,100,-0.2,0.2);
        hi_vxy_neg.setTitleX("Vx (cm)");
        hi_vxy_neg.setTitleY("Vy (cm)"); 
        H2F hi_vz_vs_theta_pos = new H2F("hi_vz_vs_theta_pos","hi_vz_vs_theta_pos",100, 30.,120.,100,-15.,15);
        hi_vz_vs_theta_pos.setTitleX("#theta (deg)");
        hi_vz_vs_theta_pos.setTitleY("Vz (cm)");
        H2F hi_vz_vs_theta_neg = new H2F("hi_vz_vs_theta_neg","hi_vz_vs_theta_neg",100, 30.,120.,100,-15.,15);
        hi_vz_vs_theta_neg.setTitleX("#theta (deg)");
        hi_vz_vs_theta_neg.setTitleY("Vz (cm)");
        H2F hi_vz_vs_phi_pos = new H2F("hi_vz_vs_phi_pos","hi_vz_vs_phi_pos",200,-15.,15.,200,-180,180);
        hi_vz_vs_phi_pos.setTitleX("Vz (cm)");
        hi_vz_vs_phi_pos.setTitleY("#phi (deg)");
        H2F hi_vz_vs_phi_neg = new H2F("hi_vz_vs_phi_neg","hi_vz_vs_phi_neg",200,-5.,5.,200,-180,180);
        hi_vz_vs_phi_neg.setTitleX("Vz (cm)");
        hi_vz_vs_phi_neg.setTitleY("#phi (deg)");
        DataGroup vertex = new DataGroup(3,2);
        vertex.addDataSet(hi_vz_vs_theta_pos, 0);
        vertex.addDataSet(hi_vxy_pos, 1);
        vertex.addDataSet(hi_vz_vs_phi_pos, 2);
        vertex.addDataSet(hi_vz_vs_theta_neg, 3);
        vertex.addDataSet(hi_vxy_neg, 4);
        vertex.addDataSet(hi_vz_vs_phi_neg, 5);
        this.getDataGroup().add(vertex, 5);
        // trajectories
        H1F hi_p_pos_traj = new H1F("hi_p_pos_traj", "hi_p_pos_traj", 100, 0.0, 3.0);     
        hi_p_pos_traj.setTitleX("p (GeV)");
        hi_p_pos_traj.setTitleY("Counts");
        H1F hi_theta_pos_traj = new H1F("hi_theta_pos_traj", "hi_theta_pos_traj", 100, 30.0, 150.0); 
        hi_theta_pos_traj.setTitleX("#theta (deg)");
        hi_theta_pos_traj.setTitleY("Counts");
        H1F hi_path_pos_traj = new H1F("hi_path_pos_traj", "hi_path_pos_traj", 100, 30.0, 100.0); 
        hi_path_pos_traj.setTitleX("Path (cm)");
        hi_path_pos_traj.setTitleY("Counts");
        H1F hi_dpath_pos_traj = new H1F("hi_dpath_pos_traj", "hi_dpath_pos_traj", 100, 0.0, 20.0); 
        hi_dpath_pos_traj.setTitleX("#DeltaPath (cm)");
        hi_dpath_pos_traj.setTitleY("Counts");
        H1F hi_p_neg_traj = new H1F("hi_p_neg_traj", "hi_p_neg_traj", 100, 0.0, 3.0);     
        hi_p_neg_traj.setTitleX("p (GeV)");
        hi_p_neg_traj.setTitleY("Counts");
        H1F hi_theta_neg_traj = new H1F("hi_theta_neg_traj", "hi_theta_neg_traj", 100, 30.0, 150.0); 
        hi_theta_neg_traj.setTitleX("#theta (deg)");
        hi_theta_neg_traj.setTitleY("Counts");
        H1F hi_path_neg_traj = new H1F("hi_path_neg_traj", "hi_path_neg_traj", 100, 30.0, 100.0); 
        hi_path_neg_traj.setTitleX("Path (cm)");
        hi_path_neg_traj.setTitleY("Counts");
        H1F hi_dpath_neg_traj = new H1F("hi_dpath_neg_traj", "hi_dpath_neg_traj", 100, 0.0, 20.0); 
        hi_dpath_neg_traj.setTitleX("#DeltaPath (cm)");
        hi_dpath_neg_traj.setTitleY("Counts");
        DataGroup trajectories = new DataGroup(4,2);
        trajectories.addDataSet(hi_p_pos_traj,     0);
        trajectories.addDataSet(hi_theta_pos_traj, 1);
        trajectories.addDataSet(hi_path_pos_traj,  2);
        trajectories.addDataSet(hi_dpath_pos_traj, 3);
        trajectories.addDataSet(hi_p_neg_traj,     4);
        trajectories.addDataSet(hi_theta_neg_traj, 5);
        trajectories.addDataSet(hi_path_neg_traj,  6);
        trajectories.addDataSet(hi_dpath_neg_traj, 7);
        this.getDataGroup().add(trajectories, 6);
    }
    
    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(0);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_p_pos_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(1);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_theta_pos_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(2);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_path_pos_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(3);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_dpath_pos_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(4);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_p_neg_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(5);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_theta_neg_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(6);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_path_neg_traj"));
        this.getAnalysisCanvas().getCanvas("Trajectories").cd(7);
        this.getAnalysisCanvas().getCanvas("Trajectories").draw(this.getDataGroup().getItem(6).getH1F("hi_dpath_neg_traj"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_p_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_chi2_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(6);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(7);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_ndf_neg"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_p_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_mc_p_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_mc_theta_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_phi_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_mc_phi_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_chi2_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_mc_vz_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_p_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(6);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_phi_pos"));        
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(7);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_ndf_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(0);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos"));
//        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").cd(1);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vz_vs_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(2);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vxy_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(3);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vz_vs_phi_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(4);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg"));
//        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").cd(5);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vz_vs_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(6);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(6).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vxy_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(7);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(5).getH2F("hi_vz_vs_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(0);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_p_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(1);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_phi_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(3);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_z_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(4);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_p_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(5);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(6);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").cd(7);
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").draw(this.getDataGroup().getItem(4).getH2F("hi_dp_z_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(0);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dp_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dp_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(1);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dtheta_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dtheta_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dphi_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dphi_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(3);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dvz_pos"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dvz_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(4);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dp_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dp_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(5);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dtheta_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(6);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dphi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dphi_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(7);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dvz_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getF1D("f1_dvz_neg"),"same");

        
        this.getAnalysisCanvas().getCanvas("Trajectories").update();
        this.getAnalysisCanvas().getCanvas("Negative Tracks").update();
        this.getAnalysisCanvas().getCanvas("Positive Tracks").update();
        this.getAnalysisCanvas().getCanvas("Vertex").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo 2D").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo").update();
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        Particle partGenNeg = null;
        Particle partGenPos = null;
        Particle partRecNeg = null;
        Particle partRecPos = null;
        if(event.hasBank("CVTRec::Tracks")){
            DataBank  bank = event.getBank("CVTRec::Tracks");
            DataBank  traj = null;
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
            	
                int pidCode=0;
                if(bank.getByte("q", loop)==-1) pidCode = -211;
                else if(bank.getByte("q", loop)==1) pidCode = 211;
                
                else pidCode = 22;
                Particle recParticle = new Particle(
                                          pidCode,
                                          bank.getFloat("pt", loop)*Math.cos(bank.getFloat("phi0", loop)),
                                          bank.getFloat("pt", loop)*Math.sin(bank.getFloat("phi0", loop)),
                                          bank.getFloat("pt", loop)*bank.getFloat("tandip", loop),
                                          -bank.getFloat("d0", loop)*Math.sin(bank.getFloat("phi0", loop)),
                                          bank.getFloat("d0", loop)*Math.cos(bank.getFloat("phi0", loop)),
                                          bank.getFloat("z0", loop));
                recParticle.setProperty("pt", bank.getFloat("pt", loop));
                recParticle.setProperty("chi2", bank.getFloat("chi2", loop));
                recParticle.setProperty("ndf", bank.getShort("ndf", loop));
                double[] path = new double[10];
                double trajZ = 0;
                double svtZ = 0;
                if(event.hasBank("CVTRec::Trajectory")) {
                    traj = event.getBank("CVTRec::Trajectory");
                    for(int j=0; j<traj.rows(); j++) {
                        if(bank.getShort("ID",loop)==traj.getShort("id",j)) {
                            int detector = traj.getByte("detector",j);
                            int layer    = traj.getByte("layer",j);
                            if(detector==DetectorType.CVT.getDetectorId() && layer==1) {
                                svtZ = traj.getFloat("z", j);
                            }
                            if(detector==DetectorType.CVT.getDetectorId() && layer>6) {
                                path[-7+layer]= traj.getFloat("path", j);
                            }
                            else if(detector==DetectorType.CTOF.getDetectorId()) {
                                path[5+layer]= traj.getFloat("path", j);
                            }
                            else if(detector==DetectorType.CND.getDetectorId()) {
                                path[6+layer]= traj.getFloat("path", j);
                                if(layer==3) trajZ = traj.getFloat("z", j);
                            }
                        }
                    }
                }
                if(recParticle.charge()>0 && recParticle.getProperty("ndf")>=-3 && recParticle.getProperty("pt")>0.2) {
                    this.getDataGroup().getItem(2).getH1F("hi_p_pos").fill(recParticle.p());
                    this.getDataGroup().getItem(2).getH1F("hi_theta_pos").fill((float) Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_phi_pos").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(2).getH1F("hi_chi2_pos").fill(recParticle.getProperty("chi2"));
                    this.getDataGroup().getItem(2).getH1F("hi_vz_pos").fill(recParticle.vz());
                    this.getDataGroup().getItem(2).getH2F("hi_theta_p_pos").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH2F("hi_theta_phi_pos").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_ndf_pos").fill(recParticle.getProperty("ndf"));
                    this.getDataGroup().getItem(5).getH2F("hi_vz_vs_theta_pos").fill(Math.toDegrees(recParticle.theta()),recParticle.vz());
                    this.getDataGroup().getItem(5).getH2F("hi_vz_vs_phi_pos").fill(recParticle.vz(),Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(5).getH2F("hi_vxy_pos").fill(recParticle.vx(),recParticle.vy());
                    if(path[9]>0) {
                       this.getDataGroup().getItem(6).getH1F("hi_p_pos_traj").fill(recParticle.p());
                       this.getDataGroup().getItem(6).getH1F("hi_theta_pos_traj").fill(Math.toDegrees(recParticle.theta()));
                       this.getDataGroup().getItem(6).getH1F("hi_path_pos_traj").fill(path[9]);
                       for(int k=0; k<9; k++) this.getDataGroup().getItem(6).getH1F("hi_dpath_pos_traj").fill(path[k+1]-path[k]);
                    }
                }
                else if(recParticle.charge()<0 && recParticle.getProperty("ndf")>=-3 && recParticle.getProperty("pt")>0.2){
                    this.getDataGroup().getItem(1).getH1F("hi_p_neg").fill(recParticle.p());
                    this.getDataGroup().getItem(1).getH1F("hi_theta_neg").fill(Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH1F("hi_phi_neg").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(1).getH1F("hi_chi2_neg").fill(recParticle.getProperty("chi2"));
                    this.getDataGroup().getItem(1).getH1F("hi_vz_neg").fill(recParticle.vz());
                    this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));                    
                    this.getDataGroup().getItem(1).getH1F("hi_ndf_neg").fill(recParticle.getProperty("ndf"));
                    this.getDataGroup().getItem(5).getH2F("hi_vz_vs_theta_neg").fill(Math.toDegrees(recParticle.theta()),recParticle.vz());
                    this.getDataGroup().getItem(5).getH2F("hi_vz_vs_phi_neg").fill(recParticle.vz(),Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(5).getH2F("hi_vxy_neg").fill(recParticle.vx(),recParticle.vy());
                    if(path[9]>0) {
                       this.getDataGroup().getItem(6).getH1F("hi_p_neg_traj").fill(recParticle.p());
                       this.getDataGroup().getItem(6).getH1F("hi_theta_neg_traj").fill(Math.toDegrees(recParticle.theta()));
                       this.getDataGroup().getItem(6).getH1F("hi_path_neg_traj").fill(path[9]);
                       for(int k=0; k<9; k++) this.getDataGroup().getItem(6).getH1F("hi_dpath_neg_traj").fill(path[k+1]-path[k]);
                    }
                }
                if(partRecNeg==null && recParticle.charge()<0) partRecNeg=recParticle;
                if(partRecPos==null && recParticle.charge()>0) partRecPos=recParticle;
                
//                System.out.println(" q = "+recParticle.charge());
            }
        }
        if(event.hasBank("MC::Particle")) {
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
                if(genPart.charge()>0) {
                    this.getDataGroup().getItem(2).getH1F("hi_mc_p_pos").fill(genPart.p());
                    this.getDataGroup().getItem(2).getH1F("hi_mc_theta_pos").fill((float) Math.toDegrees(genPart.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_mc_phi_pos").fill(Math.toDegrees(genPart.phi()));
                    this.getDataGroup().getItem(2).getH1F("hi_mc_vz_pos").fill(genPart.vz());
                }
            }
        }
        if(event.hasBank("CVTRec::Tracks")){
            DataBank genBank = null;
            if(event.hasBank("MC::Particle")) {
                genBank = event.getBank("MC::Particle");           
    //            System.out.println(" Getting Truth info....");             
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
                    //
                    double svtZ = 0;
                    if(event.hasBank("CVTRec::Trajectory")) {
                        DataBank traj = event.getBank("CVTRec::Trajectory");
                        for(int j=0; j<traj.rows(); j++) {
                            if(event.getBank("CVTRec::Tracks").getShort("ID",loop)==traj.getShort("id",j)) {
                                int detector = traj.getByte("detector",j);
                                int layer    = traj.getByte("layer",j);
                                if(detector==DetectorType.CVT.getDetectorId() && layer==1) {
                                    svtZ = traj.getFloat("z", j);
                                }
                            }
                        }
                    }
                    //if(genPart.pid()==11  && partGenNeg==null) partGenNeg=genPart;
                    //if(genPart.pid()!=211 && partGenPos==null) partGenPos=genPart;
//                    partGenPos=genPart;
//                    System.out.println(" --> q+ = "+partRecPos.charge()+" "+((partRecPos != null)) +" "+partRecPos.p() );
//                    System.out.println(" --> q- = "+partRecNeg.charge()+" "+((partRecNeg != null)) );
                    if(partRecNeg != null && genPart.charge()==-1 && partRecNeg.getProperty("ndf")>=-3 && partRecNeg.getProperty("pt")>0.2) { 
                        //System.out.println(" Filling negative track info ...");
                        this.getDataGroup().getItem(3).getH1F("hi_dp_neg").fill((partRecNeg.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg").fill(Math.toDegrees(partRecNeg.theta()-genPart.theta()));
                        this.getDataGroup().getItem(3).getH1F("hi_dphi_neg").fill(Math.toDegrees(partRecNeg.phi()-genPart.phi()));
                        this.getDataGroup().getItem(3).getH1F("hi_dvz_neg").fill(partRecNeg.vz()-genPart.vz());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_p_neg").fill(genPart.p(),(partRecNeg.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_theta_neg").fill(Math.toDegrees(genPart.theta()),(partRecNeg.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_phi_neg").fill(Math.toDegrees(genPart.phi()),(partRecNeg.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_z_neg").fill(genPart.vz(),(partRecNeg.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_z_pos").fill(svtZ,(partRecNeg.p()-genPart.p())/genPart.p());
                    }
                    if(partRecPos != null && genPart.charge()==1 && partRecPos.getProperty("ndf")>=-3 && partRecPos.getProperty("pt")>0.2) { 
                        //System.out.println(" Filling positive track info ...");
                        this.getDataGroup().getItem(3).getH1F("hi_dp_pos").fill((partRecPos.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(3).getH1F("hi_dtheta_pos").fill(Math.toDegrees(partRecPos.theta()-genPart.theta()));
                        this.getDataGroup().getItem(3).getH1F("hi_dphi_pos").fill(Math.toDegrees(partRecPos.phi()-genPart.phi()));
                        this.getDataGroup().getItem(3).getH1F("hi_dvz_pos").fill(partRecPos.vz()-genPart.vz());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_p_pos").fill(genPart.p(),(partRecPos.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_theta_pos").fill(Math.toDegrees(genPart.theta()),(partRecPos.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_phi_pos").fill(Math.toDegrees(genPart.phi()),(partRecPos.p()-genPart.p())/genPart.p());
                        this.getDataGroup().getItem(4).getH2F("hi_dp_z_pos").fill(genPart.vz(),(partRecPos.p()-genPart.p())/genPart.p());
                    }
                }
            }
        }
    }


    @Override
    public void timerUpdate() {
//        System.out.println("Updating HBT");
        // fitting MC comparisons
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dp_pos"),     this.getDataGroup().getItem(3).getF1D("f1_dp_pos"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dtheta_pos"), this.getDataGroup().getItem(3).getF1D("f1_dtheta_pos"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dphi_pos"),   this.getDataGroup().getItem(3).getF1D("f1_dphi_pos"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dvz_pos"),    this.getDataGroup().getItem(3).getF1D("f1_dvz_pos"));     
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dp_neg"),     this.getDataGroup().getItem(3).getF1D("f1_dp_neg"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg"), this.getDataGroup().getItem(3).getF1D("f1_dtheta_neg"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dphi_neg"),   this.getDataGroup().getItem(3).getF1D("f1_dphi_neg"));
        this.fitMC(this.getDataGroup().getItem(3).getH1F("hi_dvz_neg"),    this.getDataGroup().getItem(3).getF1D("f1_dvz_neg"));     
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

}
