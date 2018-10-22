/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.pdg.PhysicsConstants;
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
public class EBHBmonitor extends AnalysisMonitor {
    
    double rfPeriod = 4.008;
    int    nFT = 0;
    int    nFTel = 0;
    
    public EBHBmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Start time", "Electrons", "Pi0", "Vertex time", "Mass", "Beta");
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
   
        // Start time
        DataGroup dc_start = new DataGroup(2,2);
        H1F hi_start_all = new H1F("hi_start_all", "hi_start_all", 300, 0.,300.);  
        hi_start_all.setTitleX("Start time (ns)"); 
        hi_start_all.setTitleY("Counts");
        dc_start.addDataSet(hi_start_all, 0);
        H1F hi_start_el = new H1F("hi_start_el", "hi_start_el", 300, 0.,300.);  
        hi_start_el.setTitleX("Start time (ns)"); 
        hi_start_el.setTitleY("Counts");
        hi_start_el.setLineColor(2);
        H1F hi_start_oth = new H1F("hi_start_oth", "hi_start_oth", 300, 0.,300.);  
        hi_start_oth.setTitleX("Start time (ns)"); 
        hi_start_oth.setTitleY("Counts");
        hi_start_oth.setLineColor(4);
        H1F hi_start_opp = new H1F("hi_start_opp", "hi_start_opp", 300, 0.,300.);  
        hi_start_opp.setTitleX("Start time (ns)"); 
        hi_start_opp.setTitleY("Counts");
        hi_start_opp.setLineColor(3);
        H1F hi_vt_el = new H1F("hi_vt_el", "hi_vt_el", 300, -rfPeriod/2, rfPeriod/2);  
        hi_vt_el.setTitleX("Electron vertex time (ns)"); 
        hi_vt_el.setTitleY("Counts");
        hi_vt_el.setFillColor(4);
        F1D f1_el = new F1D("f1_el", "[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1_el.setParameter(0, 0);
        f1_el.setParameter(1, 0);
        f1_el.setParameter(2, 0.2);
        f1_el.setLineWidth(2);
        f1_el.setLineColor(2);
        f1_el.setOptStat("1111");
        H1F hi_vt_pi = new H1F("hi_vt_pi", "hi_vt_pi", 300, -rfPeriod/2, rfPeriod/2);  
        hi_vt_pi.setTitleX("Pion vertex time (ns)"); 
        hi_vt_pi.setTitleY("Counts");
        hi_vt_pi.setFillColor(3);
        F1D f1_pi = new F1D("f1_pi", "[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1_pi.setParameter(0, 0);
        f1_pi.setParameter(1, 0);
        f1_pi.setParameter(2, 0.2);
        f1_pi.setLineWidth(2);
        f1_pi.setLineColor(2);
        f1_pi.setOptStat("1111");
        H1F hi_vt_pr = new H1F("hi_vt_pr", "hi_vt_pr", 300, -rfPeriod/2, rfPeriod/2);  
        hi_vt_pr.setTitleX("Proton vertex time (ns)"); 
        hi_vt_pr.setTitleY("Counts");
        hi_vt_pr.setFillColor(2);
        F1D f1_pr = new F1D("f1_pr", "[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
        f1_pr.setParameter(0, 0);
        f1_pr.setParameter(1, 0);
        f1_pr.setParameter(2, 0.2);
        f1_pr.setLineWidth(2);
        f1_pr.setLineColor(2);
        f1_pr.setOptStat("1111");        dc_start.addDataSet(hi_start_all, 0);
        dc_start.addDataSet(hi_start_el,  0);
        dc_start.addDataSet(hi_start_oth, 0);
        dc_start.addDataSet(hi_start_opp, 0);
        dc_start.addDataSet(hi_vt_el,     1);
        dc_start.addDataSet(f1_el,        1);
        dc_start.addDataSet(hi_vt_pi,     2);
        dc_start.addDataSet(f1_pi,        2);
        dc_start.addDataSet(hi_vt_pr,     3);
        dc_start.addDataSet(f1_pr,        3);
        this.getDataGroup().add(dc_start, 0);  
        // Electrons
        DataGroup dc_electron = new DataGroup(2,2);
        H2F hi_el_p_theta_fd = new H2F("hi_el_p_theta_fd", "hi_el_p_theta_fd", 100, 0.,10., 100, 5, 40);  
        hi_el_p_theta_fd.setTitleX("p(GeV)"); 
        hi_el_p_theta_fd.setTitleY("#theta (deg)");
        H2F hi_el_p_theta_ft = new H2F("hi_el_p_theta_ft", "hi_el_p_theta_ft", 100, 0.,10., 100, 2, 5);  
        hi_el_p_theta_ft.setTitleX("p(GeV)"); 
        hi_el_p_theta_ft.setTitleY("#theta (deg)");
        H1F hi_el_p_fd = new H1F("hi_el_p_fd", "hi_el_p_fd", 200, 0.,10.);  
        hi_el_p_fd.setTitleX("p(GeV)"); 
        hi_el_p_fd.setTitleY("Counts");
        hi_el_p_fd.setFillColor(4);
        H1F hi_el_p_ft = new H1F("hi_el_p_ft", "hi_el_p_ft", 200, 0.,10.);  
        hi_el_p_ft.setTitleX("p(GeV)"); 
        hi_el_p_ft.setTitleY("Counts");
        hi_el_p_ft.setFillColor(4);
        dc_electron.addDataSet(hi_el_p_theta_fd, 0);
        dc_electron.addDataSet(hi_el_p_theta_ft, 1);
        dc_electron.addDataSet(hi_el_p_fd, 2);
        dc_electron.addDataSet(hi_el_p_ft, 3);
        this.getDataGroup().add(dc_electron, 1);          
        // Electrons
        DataGroup dc_pi0 = new DataGroup(2,2);
        H2F hi_pi0_angle_fd = new H2F("hi_pi0_angle_fd", "hi_pi0_angle_fd",100, 0., 30., 100, 0, 3);   
        hi_pi0_angle_fd.setTitleX("#theta (deg)"); 
        hi_pi0_angle_fd.setTitleY("E(GeV)");
        H2F hi_pi0_angle_ft = new H2F("hi_pi0_angle_ft", "hi_pi0_angle_ft", 100, 0., 300., 100, 0., 6.); 
        hi_pi0_angle_ft.setTitleX("Mass (MeV)"); 
        hi_pi0_angle_ft.setTitleY("#theta (deg)");
        H1F hi_pi0_mass_fd = new H1F("hi_pi0_mass_fd", "hi_pi0_mass_fd", 200, 0.,300.);  
        hi_pi0_mass_fd.setTitleX("p(GeV)"); 
        hi_pi0_mass_fd.setTitleY("Counts");
        hi_pi0_mass_fd.setFillColor(3);
        H1F hi_pi0_mass_ft = new H1F("hi_pi0_mass_ft", "hi_pi0_mass_ft", 200, 0.,300.);  
        hi_pi0_mass_ft.setTitleX("p(GeV)"); 
        hi_pi0_mass_ft.setTitleY("Counts");
        hi_pi0_mass_ft.setFillColor(3);
        dc_pi0.addDataSet(hi_pi0_angle_fd, 0);
        dc_pi0.addDataSet(hi_pi0_angle_ft, 1);
        dc_pi0.addDataSet(hi_pi0_mass_fd, 2);
        dc_pi0.addDataSet(hi_pi0_mass_ft, 3);
        this.getDataGroup().add(dc_pi0, 2);          
        // Vertex time
        DataGroup dc_time = new DataGroup(2,2);
        H2F hi_time_pos_ftof = new H2F("hi_time_pos_ftof", "hi_time_pos_ftof", 100, 0., 5., 100, -5, 5);  
        hi_time_pos_ftof.setTitleX("p (GeV)"); 
        hi_time_pos_ftof.setTitleY("Vertex time (ns)");
        H2F hi_time_neg_ftof = new H2F("hi_time_neg_ftof", "hi_time_neg_ftof", 100, 0., 5., 100, -5, 5); 
        hi_time_neg_ftof.setTitleX("p (GeV)"); 
        hi_time_neg_ftof.setTitleY("Vertex time (ns)");
        H2F hi_time_pos_ctof = new H2F("hi_time_pos_ctof", "hi_time_pos_ctof", 100, 0., 3., 100, -30, 30); 
        hi_time_pos_ctof.setTitleX("p (GeV)"); 
        hi_time_pos_ctof.setTitleY("Vertex time (ns)");
        H2F hi_time_neg_ctof = new H2F("hi_time_neg_ctof", "hi_time_neg_ctof", 100, 0., 3., 100, -30, 30);  
        hi_time_neg_ctof.setTitleX("p (GeV)"); 
        hi_time_neg_ctof.setTitleY("Vertex time (ns)");
        dc_time.addDataSet(hi_time_pos_ftof,  0);
        dc_time.addDataSet(hi_time_neg_ftof,  1);
        dc_time.addDataSet(hi_time_pos_ctof,  2);
        dc_time.addDataSet(hi_time_neg_ctof,  3);
        this.getDataGroup().add(dc_time, 3);  
        // mass
        DataGroup dc_mass = new DataGroup(2,2);
        H1F hi_mass_pos_ftof = new H1F("hi_mass_pos_ftof", "hi_mass_pos_ftof", 100, -1., 3.);  
        hi_mass_pos_ftof.setTitleX("Mass^2 (GeV)"); 
        hi_mass_pos_ftof.setTitleY("Counts");
        hi_mass_pos_ftof.setFillColor(32);
        H1F hi_mass_neg_ftof = new H1F("hi_mass_neg_ftof", "hi_mass_neg_ftof", 100, -1., 3.);  
        hi_mass_neg_ftof.setTitleX("Mass^2 (GeV)"); 
        hi_mass_neg_ftof.setTitleY("Counts");
        hi_mass_neg_ftof.setFillColor(34);
        H1F hi_mass_pos_ctof = new H1F("hi_mass_pos_ctof", "hi_mass_pos_ctof", 100, -1., 3.);  
        hi_mass_pos_ctof.setTitleX("Mass^2 (GeV)"); 
        hi_mass_pos_ctof.setTitleY("Counts");
        hi_mass_pos_ctof.setFillColor(32);
        H1F hi_mass_neg_ctof = new H1F("hi_mass_neg_ctof", "hi_mass_neg_ctof", 100, -1., 3.);  
        hi_mass_neg_ctof.setTitleX("Mass^2 (GeV)"); 
        hi_mass_neg_ctof.setTitleY("Counts");
        hi_mass_neg_ctof.setFillColor(34);
        dc_mass.addDataSet(hi_mass_pos_ftof, 0);
        dc_mass.addDataSet(hi_mass_neg_ftof, 1);
        dc_mass.addDataSet(hi_mass_pos_ctof, 2);
        dc_mass.addDataSet(hi_mass_neg_ctof, 3);
        this.getDataGroup().add(dc_mass, 4);          
        // beta
        DataGroup dc_beta = new DataGroup(2,2);
        H2F hi_beta_pos_ftof = new H2F("hi_beta_pos_ftof", "hi_beta_pos_ftof", 100, 0., 5., 100, 0., 1.5);  
        hi_beta_pos_ftof.setTitleX("p (GeV)"); 
        hi_beta_pos_ftof.setTitleY("#beta");
        H2F hi_beta_neg_ftof = new H2F("hi_beta_neg_ftof", "hi_beta_neg_ftof", 100, 0., 5., 100, 0., 1.5);  
        hi_beta_neg_ftof.setTitleX("p (GeV)"); 
        hi_beta_neg_ftof.setTitleY("#beta");
        H2F hi_beta_pos_ctof = new H2F("hi_beta_pos_ctof", "hi_beta_pos_ctof", 100, 0., 3., 100, 0.2, 1.5);  
        hi_beta_pos_ctof.setTitleX("p (GeV)"); 
        hi_beta_pos_ctof.setTitleY("#beta");
        H2F hi_beta_neg_ctof = new H2F("hi_beta_neg_ctof", "hi_beta_neg_ctof", 100, 0., 3., 100, 0.2, 1.5);  
        hi_beta_neg_ctof.setTitleX("p (GeV)"); 
        hi_beta_neg_ctof.setTitleY("#beta");
        F1D fpion = new F1D("fpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 5.0);
        F1D fkaon = new F1D("fkaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 5.0);
        F1D fprot = new F1D("fprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 5.0);
        F1D cpion = new F1D("cpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 3.0);
        F1D ckaon = new F1D("ckaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 3.0);
        F1D cprot = new F1D("cprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 3.0);
        dc_beta.addDataSet(hi_beta_pos_ftof,  0);
        dc_beta.addDataSet(fpion,  0);
        dc_beta.addDataSet(fkaon,  0);
        dc_beta.addDataSet(fprot,  0);
        dc_beta.addDataSet(hi_beta_neg_ftof,  1);
        dc_beta.addDataSet(fpion,  1);
        dc_beta.addDataSet(fkaon,  1);
        dc_beta.addDataSet(fprot,  1);
        dc_beta.addDataSet(hi_beta_pos_ctof,  2);
        dc_beta.addDataSet(cpion,  2);
        dc_beta.addDataSet(ckaon,  2);
        dc_beta.addDataSet(cprot,  2);
        dc_beta.addDataSet(hi_beta_neg_ctof,  3);
        dc_beta.addDataSet(cpion,  3);
        dc_beta.addDataSet(ckaon,  3);
        dc_beta.addDataSet(cprot,  3);        
        this.getDataGroup().add(dc_beta, 5);  
    }
    
    @Override
    public void plotHistos() {

        this.getAnalysisCanvas().getCanvas("Beta").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Beta").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Beta").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Mass").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Mass").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Mass").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Vertex time").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Vertex time").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Vertex time").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Pi0").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Pi0").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Pi0").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Electrons").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Start time").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Start time").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Start time").setGridY(false);
        
        
        this.getAnalysisCanvas().getCanvas("Beta").cd(0);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fkaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(1);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_neg_ftof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fkaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("fprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(2);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(3);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_neg_ctof"));    
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(5).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Mass").cd(0);
        this.getAnalysisCanvas().getCanvas("Mass").getPad(0).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Mass").cd(1);
        this.getAnalysisCanvas().getCanvas("Mass").getPad(1).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_neg_ftof"));
        this.getAnalysisCanvas().getCanvas("Mass").cd(2);
        this.getAnalysisCanvas().getCanvas("Mass").getPad(2).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Mass").cd(3);
        this.getAnalysisCanvas().getCanvas("Mass").getPad(3).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_neg_ctof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(0);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(3).getH2F("hi_time_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(1);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(3).getH2F("hi_time_neg_ftof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(2);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(3).getH2F("hi_time_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(3);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(3).getH2F("hi_time_neg_ctof"));
        this.getAnalysisCanvas().getCanvas("Pi0").cd(0);
        this.getAnalysisCanvas().getCanvas("Pi0").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Pi0").draw(this.getDataGroup().getItem(2).getH2F("hi_pi0_angle_fd"));
        this.getAnalysisCanvas().getCanvas("Pi0").cd(1);
        this.getAnalysisCanvas().getCanvas("Pi0").draw(this.getDataGroup().getItem(2).getH2F("hi_pi0_angle_ft"));
        this.getAnalysisCanvas().getCanvas("Pi0").cd(2);
        this.getAnalysisCanvas().getCanvas("Pi0").draw(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_fd"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(0);
        this.getAnalysisCanvas().getCanvas("Pi0").cd(3);
        this.getAnalysisCanvas().getCanvas("Pi0").draw(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_ft"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(0);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH2F("hi_el_p_theta_fd"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(1);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH2F("hi_el_p_theta_ft"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(2);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH1F("hi_el_p_fd"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(3);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH1F("hi_el_p_ft"));
        this.getAnalysisCanvas().getCanvas("Start time").cd(0);
//        this.getAnalysisCanvas().getCanvas("Start time").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_start_all"));
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_start_el"),"same");
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_start_oth"),"same");
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_start_opp"),"same");
        this.getAnalysisCanvas().getCanvas("Start time").cd(1);
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_vt_el"));
        this.getAnalysisCanvas().getCanvas("Start time").cd(2);
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_vt_pi"));
        this.getAnalysisCanvas().getCanvas("Start time").cd(3);
        this.getAnalysisCanvas().getCanvas("Start time").draw(this.getDataGroup().getItem(0).getH1F("hi_vt_pr"));
     
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recDeteEB = null;
        DataBank recFTagEB = null;
        DataBank recEvenEB = null;
        DataBank recRunRF  = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("RECHB::Particle"))        recBankEB   = event.getBank("RECHB::Particle");
        if(event.hasBank("RECHB::Scintillator"))    recDeteEB   = event.getBank("RECHB::Scintillator");
        if(event.hasBank("RECHB::ForwardTagger"))   recFTagEB   = event.getBank("RERECHBC::ForwardTagger");
        if(event.hasBank("RECHB::Event"))           recEvenEB   = event.getBank("RECHB::Event");
        if(event.hasBank("RUN::rf"))                recRunRF    = event.getBank("RUN::rf");
        int ev = 0;
        if(recRun != null) ev=recRun.getInt("event",0);
//        System.out.println(ev); 

        // Decoding Trigger Bits
        boolean[] trigger_bits = new boolean[32];
        if (event.hasBank("RUN::config")) {
            DataBank bank = event.getBank("RUN::config");
            long TriggerWord = bank.getLong("trigger", 0) & 0xFFFFFFFF;
            for (int i = 31; i >= 0; i--) {
                trigger_bits[i] = (TriggerWord & (1 << i)) != 0;
            }
        }
        if(recEvenEB!=null && recBankEB!=null) {
            ArrayList<Particle> gammasFD = new ArrayList();
            ArrayList<Particle> gammasFT = new ArrayList();
            int nrows = recBankEB.rows();
            for(int i=0; i<nrows; i++) {
                int charge = recBankEB.getByte("charge", i);
                float px = recBankEB.getFloat("px", i);
                float py = recBankEB.getFloat("py", i);
                float pz = recBankEB.getFloat("pz", i);
                short status = recBankEB.getShort("status", i);
                if(status>=1000 && status<2000) {
                    nFT++;
//                    if(charge==-1) nFTel++;
                }
                if(charge==0) {
                    Particle recParticle = new Particle(22,px,py,pz,0,0,0);
                    if(status>2000 && status<3000 && (status-2000)<100 && (status-2000)>10) {
                        gammasFD.add(recParticle);
                    }
                    else if(status>=1000 && status<2000) {
                         gammasFT.add(recParticle);                    
                    }
                }
            }
            if (gammasFT.size() >= 2) {
                for (int i1 = 0; i1 < gammasFT.size(); i1++) {
                    for (int i2 = i1 + 1; i2 < gammasFT.size(); i2++) {
                        Particle partGamma1 = gammasFT.get(i1);
                        Particle partGamma2 = gammasFT.get(i2);
                        Particle partPi0 = new Particle();
                        partPi0.copy(partGamma1);
                        partPi0.combine(partGamma2, +1);
                        double invmass = Math.sqrt(partPi0.mass2());
                        double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
                        if (angle > 1.5) {
                            this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_ft").fill(invmass * 1000);
                        }
                        this.getDataGroup().getItem(2).getH2F("hi_pi0_angle_ft").fill(invmass * 1000, angle);
                    }
                }
            }
            if (gammasFD.size() >= 2) {
                for (int i1 = 0; i1 < gammasFD.size(); i1++) {
                    for (int i2 = i1 + 1; i2 < gammasFD.size(); i2++) {
                        Particle partGamma1 = gammasFD.get(i1);
                        Particle partGamma2 = gammasFD.get(i2);
                        Particle partPi0 = new Particle();
                        partPi0.copy(partGamma1);
                        partPi0.combine(partGamma2, +1);
                        double invmass = Math.sqrt(partPi0.mass2());
                        double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
//                        if(Math.sqrt(partGamma1.p() + partGamma2.p()>0.1){
                        if(angle>3 && angle<17)  {
                            this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_fd").fill(invmass * 1000);
                        }
                        this.getDataGroup().getItem(2).getH2F("hi_pi0_angle_fd").fill(angle, Math.sqrt(partGamma1.p() * partGamma2.p()));
                    }
                }
            }
        }
        if(recEvenEB!=null && recBankEB!=null) {
            double startTime = recEvenEB.getFloat("STTime", 0);
            double    rfTime = recEvenEB.getFloat("RFTime", 0);
//            System.out.println(recEvenEB.getDouble("LT", 0) + " " + recEvenEB.getFloat("BCG", 0));
            this.getDataGroup().getItem(0).getH1F("hi_start_all").fill(startTime);
            int nrows = recBankEB.rows();
            int trigger = 0;
            for(int i=0; i<nrows; i++) {
                int pid = recBankEB.getInt("pid", i);
                float px = recBankEB.getFloat("px", i);
                float py = recBankEB.getFloat("py", i);
                float pz = recBankEB.getFloat("pz", i);
                float vx = recBankEB.getFloat("vx", i);
                float vy = recBankEB.getFloat("vy", i);
                float vz = recBankEB.getFloat("vz", i);
                float beta = recBankEB.getFloat("beta", i);
                short status = recBankEB.getShort("status", i);
                Particle recParticle = null;
                if(pid!=0) {
                    recParticle = new Particle(pid,px,py,pz,vx,vy,vz);
                    recParticle.setProperty("status", (double) status);
                }
                if(pid==11) {
                    if(recParticle.getProperty("status")>=2000) {
                        this.getDataGroup().getItem(1).getH2F("hi_el_p_theta_fd").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                        this.getDataGroup().getItem(1).getH1F("hi_el_p_fd").fill(recParticle.p());
                    }
                    else {
                        this.getDataGroup().getItem(1).getH2F("hi_el_p_theta_ft").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                        this.getDataGroup().getItem(1).getH1F("hi_el_p_ft").fill(recParticle.p());   
                        nFTel++;
                    }
                }
                if(i==0) trigger = pid;
                if(startTime>-100 && trigger==11 && i>0 && recParticle!=null) {
                    double mass2   = Math.pow(recParticle.p()/beta, 2)-recParticle.p()*recParticle.p();
                    if(recParticle.getProperty("status")>=4000) {
                        if(recParticle.charge()==1)  {
                            this.getDataGroup().getItem(5).getH2F("hi_beta_pos_ctof").fill(recParticle.p(),beta);                    
                            this.getDataGroup().getItem(4).getH1F("hi_mass_pos_ctof").fill(mass2);                    
                        }
                        if(recParticle.charge()==-1) {
                            this.getDataGroup().getItem(5).getH2F("hi_beta_neg_ctof").fill(recParticle.p(),beta);
                            this.getDataGroup().getItem(4).getH1F("hi_mass_neg_ctof").fill(mass2);                    
                        }
                    }
                    else if(recParticle.getProperty("status")>=2000) {
                        if(recParticle.charge()==1)  {
                            this.getDataGroup().getItem(5).getH2F("hi_beta_pos_ftof").fill(recParticle.p(),beta);                    
                            this.getDataGroup().getItem(4).getH1F("hi_mass_pos_ftof").fill(mass2);                    
                        }
                        if(recParticle.charge()==-1) {
                            this.getDataGroup().getItem(5).getH2F("hi_beta_neg_ftof").fill(recParticle.p(),beta);
                            this.getDataGroup().getItem(4).getH1F("hi_mass_neg_ftof").fill(mass2);                    
                        }
                    }
                }
            }
            if(trigger==11) this.getDataGroup().getItem(0).getH1F("hi_start_el").fill(startTime);
            else if(trigger!=0) this.getDataGroup().getItem(0).getH1F("hi_start_oth").fill(startTime);
            if(trigger_bits[19] || trigger_bits[20] || trigger_bits[21]) this.getDataGroup().getItem(0).getH1F("hi_start_opp").fill(startTime);
            if(recDeteEB!=null) {
                nrows = recDeteEB.rows();
                for(int i=0; i<recDeteEB.rows(); i++) {
                    int pindex   = recDeteEB.getShort("pindex", i);
                    int detector = recDeteEB.getByte("detector", i);
                    float time   = recDeteEB.getFloat("time", i);
                    float path   = recDeteEB.getFloat("path", i);
                    int pid      = recBankEB.getInt("pid", pindex);
                    int charge   = recBankEB.getByte("charge", pindex);
                    float px     = recBankEB.getFloat("px", pindex);
                    float py     = recBankEB.getFloat("py", pindex);
                    float pz     = recBankEB.getFloat("pz", pindex);
                    float vx     = recBankEB.getFloat("vx", pindex);
                    float vy     = recBankEB.getFloat("vy", pindex);
                    float vz     = recBankEB.getFloat("vz", pindex);
                    float beta   = recBankEB.getFloat("beta", pindex);
                    int status   = recBankEB.getShort("status", pindex);
                    if(charge!=0) {
                        Particle recParticle = new Particle(211*charge,px,py,pz,vx,vy,vz);
                        recParticle.setProperty("status", (double) status);
                        double betap = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+recParticle.mass2());
                        double dt   = (time - path/(betap*PhysicsConstants.speedOfLight()) - startTime);                
                        if(startTime>-100 && trigger==11 && pindex>0) {                    
                            if(recParticle.getProperty("status")>=4000) {
                                if(recParticle.charge()==1)  {
                                    this.getDataGroup().getItem(3).getH2F("hi_time_pos_ctof").fill(recParticle.p(),dt);
                                }
                                if(recParticle.charge()==-1) {
                                    this.getDataGroup().getItem(3).getH2F("hi_time_neg_ctof").fill(recParticle.p(),dt);
                                }
                            }
                            else if(recParticle.getProperty("status")>=2000) {
                                if(recParticle.charge()==1)  {
                                    this.getDataGroup().getItem(3).getH2F("hi_time_pos_ftof").fill(recParticle.p(),dt);                    
                                }
                                if(recParticle.charge()==-1) {
                                    this.getDataGroup().getItem(3).getH2F("hi_time_neg_ftof").fill(recParticle.p(),dt);
                                }
                            }
                        }
                        if(trigger==11) {
                            if(pid==11) {
                                betap = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p());
                                dt    = (time - path/(betap*PhysicsConstants.speedOfLight()) - rfTime);
                                dt    = (dt +1000.5*this.rfPeriod)%this.rfPeriod-0.5*this.rfPeriod;
                                this.getDataGroup().getItem(0).getH1F("hi_vt_el").fill(dt);
                            }
                            else if(pid==211 || pid==-211) {
                                betap = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+PhysicsConstants.massPionCharged()*PhysicsConstants.massPionCharged());
                                dt    = (time - path/(betap*PhysicsConstants.speedOfLight()) - rfTime);
                                dt    = (dt +1000.5*this.rfPeriod)%this.rfPeriod-0.5*this.rfPeriod;
                                this.getDataGroup().getItem(0).getH1F("hi_vt_pi").fill(dt);
                            }
                            else if(pid==2212) {
                                betap = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+PhysicsConstants.massProton()*PhysicsConstants.massProton());
                                dt    = (time - path/(betap*PhysicsConstants.speedOfLight()) - rfTime);
                                dt    = (dt +1000.5*this.rfPeriod)%this.rfPeriod-0.5*this.rfPeriod;
                                this.getDataGroup().getItem(0).getH1F("hi_vt_pr").fill(dt);
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

//        System.out.println(nFT + " " + nFTel);
    }
    
    @Override
    public void analyze() {
        this.fitRF(this.getDataGroup().getItem(0).getH1F("hi_vt_el"),this.getDataGroup().getItem(0).getF1D("f1_el"));
        this.fitRF(this.getDataGroup().getItem(0).getH1F("hi_vt_pi"),this.getDataGroup().getItem(0).getF1D("f1_pi"));
        this.fitRF(this.getDataGroup().getItem(0).getH1F("hi_vt_pr"),this.getDataGroup().getItem(0).getF1D("f1_pr"));
    }
    
    public void fitRF(H1F hirf,F1D f1rf) {

            double mean = hirf.getDataX(hirf.getMaximumBin());
        double amp = hirf.getBinContent(hirf.getMaximumBin());
        double sigma = hirf.getRMS();
        f1rf.setParameter(0, amp);
        f1rf.setParameter(1, mean);
        f1rf.setParameter(2, sigma);
        double rmax = Math.min(mean + 2. * Math.abs(sigma),  this.rfPeriod/2);
        double rmin = Math.max(mean - 2. * Math.abs(sigma), -this.rfPeriod/2);
        f1rf.setRange(rmin, rmax);
        DataFitter.fit(f1rf, hirf, "Q"); //No options uses error for sigma 
        hirf.setFunction(null);
        mean = f1rf.getParameter(1);
        sigma = f1rf.getParameter(2);
        rmax = Math.min(mean + 2. * Math.abs(sigma),  this.rfPeriod/2);
        rmin = Math.max(mean - 2. * Math.abs(sigma), -this.rfPeriod/2);
        f1rf.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1rf, hirf, "Q"); //No options uses error for sigma 
        hirf.setFunction(null);
    }
}
