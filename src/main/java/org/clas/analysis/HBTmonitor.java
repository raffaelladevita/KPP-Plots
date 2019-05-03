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

/**
 *
 * @author devita
 */
public class HBTmonitor extends AnalysisMonitor {
    

    public HBTmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Monte Carlo","Vertex","Positive Tracks","Negative Tracks");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // negative tracks
        H1F hi_p_neg = new H1F("hi_p_neg", "hi_p_neg", 100, 0.0, 8.0);     
        hi_p_neg.setTitleX("p (GeV)");
        hi_p_neg.setTitleY("Counts");
        H1F hi_theta_neg = new H1F("hi_theta_neg", "hi_theta_neg", 100, 0.0, 40.0); 
        hi_theta_neg.setTitleX("#theta (deg)");
        hi_theta_neg.setTitleY("Counts");
        H1F hi_phi_neg = new H1F("hi_phi_neg", "hi_phi_neg", 100, -180.0, 180.0);   
        hi_phi_neg.setTitleX("#phi (deg)");
        hi_phi_neg.setTitleY("Counts");
        H1F hi_chi2_neg = new H1F("hi_chi2_neg", "hi_chi2_neg", 100, 0.0, 180.0);   
        hi_chi2_neg.setTitleX("#chi2");
        hi_chi2_neg.setTitleY("Counts");
        H1F hi_vz_neg = new H1F("hi_vz_neg", "hi_vz_neg", 100, -15.0, 15.0);   
        hi_vz_neg.setTitleX("Vz (cm)");
        hi_vz_neg.setTitleY("Counts");
        H1F hi_vz_neg_cut = new H1F("hi_vz_neg_cut", "hi_vz_neg_cut", 100, -15.0, 15.0);   
        hi_vz_neg_cut.setTitleX("Vz (cm)");
        hi_vz_neg_cut.setTitleY("Counts");
        hi_vz_neg_cut.setLineColor(2);
        F1D f1_vz_neg = new F1D("f1_vz_neg","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1_vz_neg.setParameter(0, 0);
        f1_vz_neg.setParameter(1, 0);
        f1_vz_neg.setParameter(2, 1.0);
        f1_vz_neg.setLineWidth(2);
        f1_vz_neg.setLineColor(1);
        f1_vz_neg.setOptStat("1111");
        H2F hi_theta_p_neg = new H2F("hi_theta_p_neg", "hi_theta_p_neg", 100, 0.0, 8.0, 100, 0.0, 40.0); 
        hi_theta_p_neg.setTitleX("p (GeV)");
        hi_theta_p_neg.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_neg = new H2F("hi_theta_phi_neg", "hi_theta_phi_neg", 200, -180.0, 180.0, 200, 0.0, 40.0); 
        hi_theta_phi_neg.setTitleX("#phi (deg)");
        hi_theta_phi_neg.setTitleY("#theta (deg)");        
        H2F hi_chi2_vz_neg = new H2F("hi_chi2_vz_neg", "hi_chi2_vz_neg", 100, -15.0, 15.0, 100, 0.0, 180.0);   
        hi_chi2_vz_neg.setTitleX("Vz (cm)");
        hi_chi2_vz_neg.setTitleY("#chi2");
        DataGroup dg_neg = new DataGroup(4,2);
        dg_neg.addDataSet(hi_p_neg, 0);
        dg_neg.addDataSet(hi_theta_neg, 1);
        dg_neg.addDataSet(hi_phi_neg, 2);
        dg_neg.addDataSet(hi_chi2_neg, 3);
        dg_neg.addDataSet(hi_vz_neg, 4);
        dg_neg.addDataSet(hi_vz_neg_cut, 4);
        dg_neg.addDataSet(f1_vz_neg, 4);
        dg_neg.addDataSet(hi_theta_p_neg, 5);
        dg_neg.addDataSet(hi_theta_phi_neg, 6);
        dg_neg.addDataSet(hi_chi2_vz_neg, 7);
        this.getDataGroup().add(dg_neg, 1);
        // positive trakcs
        H1F hi_p_pos = new H1F("hi_p_pos", "hi_p_pos", 100, 0.0, 8.0);     
        hi_p_pos.setTitleX("p (GeV)");
        hi_p_pos.setTitleY("Counts");
        H1F hi_theta_pos = new H1F("hi_theta_pos", "hi_theta_pos", 100, 0.0, 40.0); 
        hi_theta_pos.setTitleX("#theta (deg)");
        hi_theta_pos.setTitleY("Counts");
        H1F hi_phi_pos = new H1F("hi_phi_pos", "hi_phi_pos", 100, -180.0, 180.0);   
        hi_phi_pos.setTitleX("#phi (deg)");
        hi_phi_pos.setTitleY("Counts");
        H1F hi_chi2_pos = new H1F("hi_chi2_pos", "hi_chi2_pos", 100, 0.0, 180.0);   
        hi_chi2_pos.setTitleX("#chi2");
        hi_chi2_pos.setTitleY("Counts");
        H1F hi_vz_pos = new H1F("hi_vz_pos", "hi_vz_pos", 100, -15.0, 15.0);   
        hi_vz_pos.setTitleX("Vz (cm)");
        hi_vz_pos.setTitleY("Counts");
        H1F hi_vz_pos_cut = new H1F("hi_vz_pos_cut", "hi_vz_pos_cut", 100, -15.0, 15.0);   
        hi_vz_pos_cut.setTitleX("Vz (cm)");
        hi_vz_pos_cut.setTitleY("Counts");
        hi_vz_pos_cut.setLineColor(2);
        F1D f1_vz_pos = new F1D("f1_vz_pos","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1_vz_pos.setParameter(0, 0);
        f1_vz_pos.setParameter(1, 0);
        f1_vz_pos.setParameter(2, 1.0);
        f1_vz_pos.setLineWidth(2);
        f1_vz_pos.setLineColor(1);
        f1_vz_pos.setOptStat("1111");
        H2F hi_theta_p_pos = new H2F("hi_theta_p_pos", "hi_theta_p_pos", 100, 0.0, 8.0, 100, 0.0, 40.0); 
        hi_theta_p_pos.setTitleX("p (GeV)");
        hi_theta_p_pos.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_pos = new H2F("hi_theta_phi_pos", "hi_theta_phi_pos", 200, -180.0, 180.0, 200, 0.0, 40.0); 
        hi_theta_phi_pos.setTitleX("#phi (deg)");
        hi_theta_phi_pos.setTitleY("#theta (deg)");  
        H2F hi_chi2_vz_pos = new H2F("hi_chi2_vz_pos", "hi_chi2_vz_pos", 100, -15.0, 15.0, 100, 0.0, 180.0);   
        hi_chi2_vz_pos.setTitleX("Vz (cm)");
        hi_chi2_vz_pos.setTitleY("#chi2");
        DataGroup dg_pos = new DataGroup(4,2);
        dg_pos.addDataSet(hi_p_pos, 0);
        dg_pos.addDataSet(hi_theta_pos, 1);
        dg_pos.addDataSet(hi_phi_pos, 2);
        dg_pos.addDataSet(hi_chi2_pos, 3);
        dg_pos.addDataSet(hi_vz_pos, 4);
        dg_pos.addDataSet(hi_vz_pos_cut, 4);
        dg_pos.addDataSet(f1_vz_pos, 4);
        dg_pos.addDataSet(hi_theta_p_pos, 5);
        dg_pos.addDataSet(hi_theta_phi_pos, 6);
        dg_pos.addDataSet(hi_chi2_vz_pos, 7);
        this.getDataGroup().add(dg_pos, 2);
        // mc comparison
        H1F hi_dp_pos = new H1F("hi_dp_pos", "hi_dp_pos", 100, -0.8, 0.8); 
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
        H1F hi_dphi_pos = new H1F("hi_dphi_pos", "hi_dphi_pos", 100, -8.0, 8.0); 
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
        H1F hi_dvz_pos = new H1F("hi_dvz_pos", "hi_dvz_pos", 100, -20.0, 20.0);   
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
        H1F hi_dp_neg = new H1F("hi_dp_neg", "hi_dp_neg", 100, -0.8, 0.8); 
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
        H1F hi_dphi_neg = new H1F("hi_dphi_neg", "hi_dphi_neg", 100, -8.0, 8.0); 
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
        H1F hi_dvz_neg = new H1F("hi_dvz_neg", "hi_dvz_neg", 100, -20.0, 20.0);   
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
        // vertex
        H2F hi_vxy_pos = new H2F("hi_vxy_pos","hi_vxy_pos",100,-15.,15.,100,-15.,15);
        hi_vxy_pos.setTitleX("Vx (cm)");
        hi_vxy_pos.setTitleY("Vy (cm)");
        H2F hi_vxy_neg = new H2F("hi_vxy_neg","hi_vxy_neg",100,-15.,15.,100,-15.,15);
        hi_vxy_neg.setTitleX("Vx (cm)");
        hi_vxy_neg.setTitleY("Vy (cm)"); 
        H2F hi_vz_vs_theta_pos = new H2F("hi_vz_vs_theta_pos","hi_vz_vs_theta_pos",100, 5.,40.,100,-15.,15);
        hi_vz_vs_theta_pos.setTitleX("#theta (deg)");
        hi_vz_vs_theta_pos.setTitleY("Vz (cm)");
        H2F hi_vz_vs_theta_neg = new H2F("hi_vz_vs_theta_neg","hi_vz_vs_theta_neg",100, 5.,40.,100,-15.,15);
        hi_vz_vs_theta_neg.setTitleX("#theta (deg)");
        hi_vz_vs_theta_neg.setTitleY("Vz (cm)");
        DataGroup vertex = new DataGroup(2,2);
        vertex.addDataSet(hi_vz_vs_theta_pos, 0);
        vertex.addDataSet(hi_vxy_pos, 1);
        vertex.addDataSet(hi_vz_vs_theta_neg, 2);
        vertex.addDataSet(hi_vxy_neg, 3);
        this.getDataGroup().add(vertex, 4);
    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Negative Tracks").divide(4,2);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").divide(4,2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").divide(4, 2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Vertex").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Vertex").setGridY(false);
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
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getF1D("f1_vz_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(6);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(7);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_chi2_vz_neg"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_p_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_phi_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_chi2_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getF1D("f1_vz_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_p_pos"));
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(6);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_phi_pos"));        
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(7);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(this.getDataGroup().getItem(2).getH2F("hi_chi2_vz_pos"));
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
        this.getAnalysisCanvas().getCanvas("Vertex").cd(0);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(2).getH1F("hi_vz_pos_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(2).getF1D("f1_vz_pos"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").cd(1);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(2);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vxy_pos"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(3);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(1).getF1D("f1_vz_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex").cd(4);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Vertex").cd(5);
        this.getAnalysisCanvas().getCanvas("Vertex").getPad(5).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex").draw(this.getDataGroup().getItem(4).getH2F("hi_vxy_neg"));
        
        this.getAnalysisCanvas().getCanvas("Negative Tracks").update();
        this.getAnalysisCanvas().getCanvas("Positive Tracks").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo").update();
        this.getAnalysisCanvas().getCanvas("Vertex").update();
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        Particle partGenNeg = null;
        Particle partGenPos = null;
        Particle partRecNeg = null;
        Particle partRecPos = null;
        if(event.hasBank("HitBasedTrkg::HBTracks")==true){
            DataBank  bank = event.getBank("HitBasedTrkg::HBTracks");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
                int pidCode=0;
                if(bank.getByte("q", loop)==-1) pidCode = 11;
                else if(bank.getByte("q", loop)==1) pidCode = 211;
                else pidCode = 22;
                Particle recParticle = new Particle(
                                          pidCode,
                                          bank.getFloat("p0_x", loop),
                                          bank.getFloat("p0_y", loop),
                                          bank.getFloat("p0_z", loop),
                                          bank.getFloat("Vtx0_x", loop),
                                          bank.getFloat("Vtx0_y", loop),
                                          bank.getFloat("Vtx0_z", loop));
                if(bank.getShort("ndf", loop)>0) recParticle.setProperty("chi2", bank.getFloat("chi2", loop)/bank.getShort("ndf", loop));
                if(recParticle.charge()>0) {
                    this.getDataGroup().getItem(2).getH1F("hi_p_pos").fill(recParticle.p());
                    this.getDataGroup().getItem(2).getH1F("hi_theta_pos").fill(Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_phi_pos").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(2).getH1F("hi_chi2_pos").fill(recParticle.getProperty("chi2"));
                    this.getDataGroup().getItem(2).getH1F("hi_vz_pos").fill(recParticle.vz());
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_pos").fill(Math.toDegrees(recParticle.theta()),recParticle.vz());
                    if(recParticle.p()>2.&& Math.toDegrees(recParticle.theta())>10.) {
                        this.getDataGroup().getItem(2).getH1F("hi_vz_pos_cut").fill(recParticle.vz());
                        this.getDataGroup().getItem(4).getH2F("hi_vxy_pos").fill(recParticle.vx(),recParticle.vy());
                    }
                    this.getDataGroup().getItem(2).getH2F("hi_theta_p_pos").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH2F("hi_theta_phi_pos").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH2F("hi_chi2_vz_pos").fill(recParticle.vz(),recParticle.getProperty("chi2"));
              }
                else {
                    this.getDataGroup().getItem(1).getH1F("hi_p_neg").fill(recParticle.p());
                    this.getDataGroup().getItem(1).getH1F("hi_theta_neg").fill(Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH1F("hi_phi_neg").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(1).getH1F("hi_chi2_neg").fill(recParticle.getProperty("chi2"));
                    this.getDataGroup().getItem(1).getH1F("hi_vz_neg").fill(recParticle.vz());
                    this.getDataGroup().getItem(4).getH2F("hi_vz_vs_theta_neg").fill(Math.toDegrees(recParticle.theta()),recParticle.vz());
                    if(recParticle.p()>2.&& Math.toDegrees(recParticle.theta())>10.) {
                        this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut").fill(recParticle.vz());
                        this.getDataGroup().getItem(4).getH2F("hi_vxy_neg").fill(recParticle.vx(),recParticle.vy());
                    }
                    this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));                    
                    this.getDataGroup().getItem(1).getH2F("hi_chi2_vz_neg").fill(recParticle.vz(),recParticle.getProperty("chi2"));
                }
                if(partRecNeg==null && recParticle.charge()<0) partRecNeg=recParticle;
                if(partRecPos==null && recParticle.charge()>0) partRecPos=recParticle;
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
                this.getDataGroup().getItem(3).getH1F("hi_dp_neg").fill((partRecNeg.p()-partGenNeg.p())/partGenNeg.p());
                this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg").fill(Math.toDegrees(partRecNeg.theta()-partGenNeg.theta()));
                this.getDataGroup().getItem(3).getH1F("hi_dphi_neg").fill(Math.toDegrees(partRecNeg.phi()-partGenNeg.phi()));
                this.getDataGroup().getItem(3).getH1F("hi_dvz_neg").fill(partRecNeg.vz()-partGenNeg.vz());
            }
            if(partGenPos != null && partRecPos != null) {
                this.getDataGroup().getItem(3).getH1F("hi_dp_pos").fill((partRecPos.p()-partGenPos.p())/partGenPos.p());
                this.getDataGroup().getItem(3).getH1F("hi_dtheta_pos").fill(Math.toDegrees(partRecPos.theta()-partGenPos.theta()));
                this.getDataGroup().getItem(3).getH1F("hi_dphi_pos").fill(Math.toDegrees(partRecPos.phi()-partGenPos.phi()));
                this.getDataGroup().getItem(3).getH1F("hi_dvz_pos").fill(partRecPos.vz()-partGenPos.vz());
            }
        }    
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating HBT");
        //fitting negative tracks vertex
        this.fitVertex(this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut"), this.getDataGroup().getItem(1).getF1D("f1_vz_neg"));
        //fitting positive tracks vertex
        this.fitVertex(this.getDataGroup().getItem(2).getH1F("hi_vz_pos_cut"), this.getDataGroup().getItem(2).getF1D("f1_vz_pos"));
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
 
    public void fitVertex(H1F hivz, F1D f1vz) {
        double mean  = hivz.getDataX(hivz.getMaximumBin());
        double amp   = hivz.getBinContent(hivz.getMaximumBin());
        double sigma = hivz.getRMS();
        if(hivz.getEntries()>500) { // first fits 
            sigma = Math.abs(f1vz.getParameter(2));       
        }
        f1vz.setParameter(0, amp);
        f1vz.setParameter(1, mean);
        f1vz.setParameter(2, sigma);
        f1vz.setRange(mean-2.*sigma,mean+2.*sigma);
        DataFitter.fit(f1vz, hivz, "Q"); //No options uses error for sigma        
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
