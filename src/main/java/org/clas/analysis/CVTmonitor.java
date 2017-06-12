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
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
/**
 *
 * @author ziegler
 */
public class CVTmonitor extends AnalysisMonitor {
    

    public CVTmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Monte Carlo","Vertex","Positive Tracks","Negative Tracks");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").divide(3,2);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").divide(3,2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").setGridY(false);
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
        H1F hi_theta_neg = new H1F("hi_theta_neg", "hi_theta_neg", 100, 50.0, 120.0); 
        hi_theta_neg.setTitleX("#theta (deg)");
        hi_theta_neg.setTitleY("Counts");
        H1F hi_phi_neg = new H1F("hi_phi_neg", "hi_phi_neg", 100, -180.0, 180.0);   
        hi_phi_neg.setTitleX("#phi (deg)");
        H1F hi_vz_neg = new H1F("hi_vz_neg", "hi_vz_neg", 100, -2.0, 2.0);   
        hi_vz_neg.setTitleX("Vz (cm)");
        hi_vz_neg.setTitleY("Counts");
        H2F hi_theta_p_neg = new H2F("hi_theta_p_neg", "hi_theta_p_neg", 100, 0.0, 3.0, 100, 50.0, 120.0); 
        hi_theta_p_neg.setTitleX("p (GeV)");
        hi_theta_p_neg.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_neg = new H2F("hi_theta_phi_neg", "hi_theta_phi_neg", 100, -180.0, 180.0, 100, 50.0, 120.0); 
        hi_theta_phi_neg.setTitleX("#phi (deg)");
        hi_theta_phi_neg.setTitleY("#theta (deg)");        
        DataGroup dg_neg = new DataGroup(3,2);
        dg_neg.addDataSet(hi_p_neg, 0);
        dg_neg.addDataSet(hi_theta_neg, 1);
        dg_neg.addDataSet(hi_phi_neg, 2);
        dg_neg.addDataSet(hi_vz_neg, 3);
        dg_neg.addDataSet(hi_theta_p_neg, 4);
        dg_neg.addDataSet(hi_theta_phi_neg, 5);
        this.getDataGroup().add(dg_neg, 1);
        // positive trakcs
        H1F hi_p_pos = new H1F("hi_p_pos", "hi_p_pos", 100, 0.0, 3.0);     
        hi_p_pos.setTitleX("p (GeV)");
        hi_p_pos.setTitleY("Counts");
        H1F hi_theta_pos = new H1F("hi_theta_pos", "hi_theta_pos", 100, 50.0, 120.0); 
        hi_theta_pos.setTitleX("#theta (deg)");
        hi_theta_pos.setTitleY("Counts");
        H1F hi_phi_pos = new H1F("hi_phi_pos", "hi_phi_pos", 100, -180.0, 180.0);   
        hi_phi_pos.setTitleX("#phi (deg)");
        H1F hi_vz_pos = new H1F("hi_vz_pos", "hi_vz_pos", 100, -2.0, 2.0);   
        hi_vz_pos.setTitleX("Vz (cm)");
        hi_vz_pos.setTitleY("Counts");
        H2F hi_theta_p_pos = new H2F("hi_theta_p_pos", "hi_theta_p_pos", 100, 0.0, 3.0, 100, 50.0, 120.0); 
        hi_theta_p_pos.setTitleX("p (GeV)");
        hi_theta_p_pos.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_pos = new H2F("hi_theta_phi_pos", "hi_theta_phi_pos", 100, -180.0, 180.0, 100, 50.0, 120.0); 
        hi_theta_phi_pos.setTitleX("#phi (deg)");
        hi_theta_phi_pos.setTitleY("#theta (deg)");  
        DataGroup dg_pos = new DataGroup(3,2);
        dg_pos.addDataSet(hi_p_pos, 0);
        dg_pos.addDataSet(hi_theta_pos, 1);
        dg_pos.addDataSet(hi_phi_pos, 2);
        dg_pos.addDataSet(hi_vz_pos, 3);
        dg_pos.addDataSet(hi_theta_p_pos, 4);
        dg_pos.addDataSet(hi_theta_phi_pos, 5);
        this.getDataGroup().add(dg_pos, 2);
        // mc comparison
        H1F hi_dp_pos = new H1F("hi_dp_pos", "hi_dp_pos", 100, -0.8, 0.8); 
        hi_dp_pos.setTitleX("#Delta P / P ");
        hi_dp_pos.setTitleY("Counts");
        hi_dp_pos.setTitle("Positive Tracks");
        H1F hi_dtheta_pos = new H1F("hi_dtheta_pos","hi_dtheta_pos", 100, -4.0, 4.0); 
        hi_dtheta_pos.setTitleX("#Delta #theta (deg)");
        hi_dtheta_pos.setTitleY("Counts");
        hi_dtheta_pos.setTitle("Positive Tracks");
        H1F hi_dphi_pos = new H1F("hi_dphi_pos", "hi_dphi_pos", 100, -8.0, 8.0); 
        hi_dphi_pos.setTitleX("#Delta #phi (deg)");
        hi_dphi_pos.setTitleY("Counts");
        hi_dphi_pos.setTitle("Positive Tracks");
        H1F hi_dvz_pos = new H1F("hi_dvz_pos", "hi_dvz_pos", 100, -20.0, 20.0);   
        hi_dvz_pos.setTitleX("#Delta Vz (cm)");
        hi_dvz_pos.setTitleY("Counts");
        hi_dvz_pos.setTitle("Positive Tracks");
        H1F hi_dp_neg = new H1F("hi_dp_neg", "hi_dp_neg", 100, -0.8, 0.8); 
        hi_dp_neg.setTitleX("#Delta P / P ");
        hi_dp_neg.setTitleY("Counts");
        hi_dp_neg.setTitle("Negative Tracks");
        H1F hi_dtheta_neg = new H1F("hi_dtheta_neg","hi_dtheta_neg", 100, -4.0, 4.0); 
        hi_dtheta_neg.setTitleX("#Delta #theta (deg)");
        hi_dtheta_neg.setTitleY("Counts");
        hi_dtheta_neg.setTitle("Negative Tracks");
        H1F hi_dphi_neg = new H1F("hi_dphi_neg", "hi_dphi_neg", 100, -8.0, 8.0); 
        hi_dphi_neg.setTitleX("#Delta #phi (deg)");
        hi_dphi_neg.setTitleY("Counts");
        hi_dphi_neg.setTitle("Negative Tracks");
        H1F hi_dvz_neg = new H1F("hi_dvz_neg", "hi_dvz_neg", 100, -20.0, 20.0);   
        hi_dvz_neg.setTitleX("#Delta Vz (cm)");
        hi_dvz_neg.setTitleY("Counts");
        hi_dvz_neg.setTitle("Negative Tracks");
        DataGroup mc = new DataGroup(4,2);
        mc.addDataSet(hi_dp_pos, 0);
        mc.addDataSet(hi_dtheta_pos, 1);
        mc.addDataSet(hi_dphi_pos, 2);
        mc.addDataSet(hi_dvz_pos, 3);
        mc.addDataSet(hi_dp_neg, 4);
        mc.addDataSet(hi_dtheta_neg, 5);
        mc.addDataSet(hi_dphi_neg, 6);
        mc.addDataSet(hi_dvz_neg, 7);
        this.getDataGroup().add(mc, 3);
        
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_p_neg);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_theta_neg);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_phi_neg);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_vz_neg);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_theta_p_neg);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Negative Tracks").draw(hi_theta_phi_neg);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_p_pos);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_theta_pos);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_phi_pos);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_vz_pos);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_theta_p_pos);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").draw(hi_theta_phi_pos);        
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(0);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dp_pos);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(1);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dtheta_pos);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dphi_pos);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(3);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dvz_pos);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(4);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dp_neg);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(5);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dtheta_neg);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(6);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dphi_neg);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(7);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(hi_dvz_neg);

        
        this.getAnalysisCanvas().getCanvas("Negative Tracks").update();
        this.getAnalysisCanvas().getCanvas("Positive Tracks").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo").update();
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        Particle partGenNeg = null;
        Particle partGenPos = null;
        Particle partRecNeg = null;
        Particle partRecPos = null;
        if(event.hasBank("CVTRec::Tracks")==true){
            DataBank  bank = event.getBank("CVTRec::Tracks");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
            	
                int pidCode=0;
                if(bank.getInt("q", loop)==-1) pidCode = 11;
                else if(bank.getInt("q", loop)==1) pidCode = 211;
                
                else pidCode = 22;
                Particle recParticle = new Particle(
                                          pidCode,
                                          bank.getFloat("pt", loop)*Math.cos(bank.getFloat("phi0", loop)),
                                          bank.getFloat("pt", loop)*Math.sin(bank.getFloat("phi0", loop)),
                                          bank.getFloat("pt", loop)*bank.getFloat("tandip", loop),
                                          -bank.getFloat("d0", loop)*Math.sin(bank.getFloat("phi0", loop)),
                                          bank.getFloat("d0", loop)*Math.cos(bank.getFloat("phi0", loop)),
                                          bank.getFloat("z0", loop));
                if(recParticle.charge()>0) {
                    this.getDataGroup().getItem(2).getH1F("hi_p_pos").fill(recParticle.p());
                    this.getDataGroup().getItem(2).getH1F("hi_theta_pos").fill((float) Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_phi_pos").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(2).getH1F("hi_vz_pos").fill(recParticle.vz());
                    this.getDataGroup().getItem(2).getH2F("hi_theta_p_pos").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(2).getH2F("hi_theta_phi_pos").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));
              }
                else {
                    this.getDataGroup().getItem(1).getH1F("hi_p_neg").fill(recParticle.p());
                    this.getDataGroup().getItem(1).getH1F("hi_theta_neg").fill(Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH1F("hi_phi_neg").fill(Math.toDegrees(recParticle.phi()));
                    this.getDataGroup().getItem(1).getH1F("hi_vz_neg").fill(recParticle.vz());
                    this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg").fill(recParticle.p(),Math.toDegrees(recParticle.theta()));
                    this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recParticle.theta()));                    
                }
                if(partRecNeg==null && recParticle.charge()<0) partRecNeg=recParticle;
                if(partRecPos==null && recParticle.charge()>0) partRecPos=recParticle;
                
//                System.out.println(" q = "+recParticle.charge());
            }
        }
        if(event.hasBank("CVTRec::Tracks")==true){
            DataBank genBank = event.getBank("MC::Particle");
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
                //if(genPart.pid()==11  && partGenNeg==null) partGenNeg=genPart;
                //if(genPart.pid()!=211 && partGenPos==null) partGenPos=genPart;
                partGenPos=genPart;
               // System.out.println(" --> q+ = "+partRecPos.charge()+" "+((partRecPos != null)) +" "+partRecPos.p() );
               // System.out.println(" --> q- = "+partRecNeg.charge()+" "+((partRecNeg != null)) );
                if(partRecNeg != null) { 
                    //System.out.println(" Filling negative track info ...");
                    this.getDataGroup().getItem(3).getH1F("hi_dp_neg").fill((partRecNeg.p()-genPart.p())/genPart.p());
                    this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg").fill(Math.toDegrees(partRecNeg.theta()-genPart.theta()));
                    this.getDataGroup().getItem(3).getH1F("hi_dphi_neg").fill(Math.toDegrees(partRecNeg.phi()-genPart.phi()));
                    this.getDataGroup().getItem(3).getH1F("hi_dvz_neg").fill(partRecNeg.vz()-genPart.vz());
                }
                if(partRecPos != null) { 
                    //System.out.println(" Filling positive track info ...");
                    this.getDataGroup().getItem(3).getH1F("hi_dp_pos").fill((partRecPos.p()-genPart.p())/genPart.p());
                    this.getDataGroup().getItem(3).getH1F("hi_dtheta_pos").fill(Math.toDegrees(partRecPos.theta()-genPart.theta()));
                    this.getDataGroup().getItem(3).getH1F("hi_dphi_pos").fill(Math.toDegrees(partRecPos.phi()-genPart.phi()));
                    this.getDataGroup().getItem(3).getH1F("hi_dvz_pos").fill(partRecPos.vz()-genPart.vz());
                }
            }
        }
    }

    @Override
    public void resetEventListener() {
        System.out.println("Resetting CVT histogram");
        this.createHistos();
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating HBT");
   }

    @Override
    public void setCanvasUpdate(int time) {
        this.getAnalysisCanvas().getCanvas("Negative Tracks").initTimer(time);
        this.getAnalysisCanvas().getCanvas("Positive Tracks").initTimer(time);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").initTimer(time);
    }
 

}
