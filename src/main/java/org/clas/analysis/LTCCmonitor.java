/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class LTCCmonitor  extends AnalysisMonitor {
        
    
    public LTCCmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Electrons");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("HTCC hits");
        summary.setFillColor(36);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        H1F hi_nphe_all = new H1F("hi_nphe_all", "hi_nphe_all", 60, 0.0, 60.0);   
        hi_nphe_all.setTitleX("N.PhE"); 
        hi_nphe_all.setTitleY("Counts"); 
        H1F hi_nphe_ele = new H1F("hi_nphe_ele", "hi_nphe_ele", 60, 0.0, 60.0);   
        hi_nphe_ele.setTitleX("N.PhE"); 
        hi_nphe_ele.setTitleY("Counts"); 
        hi_nphe_ele.setLineColor(2);
        H1F hi_time  = new H1F("hi_time", "hi_time", 200, -150., 50.);  
        hi_time.setTitleX("dT (ns)"); 
        hi_time.setTitleY("Counts");
        H2F hi_phi_TBT = new H2F("hi_phi_TBT", "hi_phi_TBT", 100, -180.0, 180.0, 100, -180.0, 180.0);  
        hi_phi_TBT.setTitleX("#phi_TBT (deg)"); 
        hi_phi_TBT.setTitleY("#phi_HTCC (deg)");
        H2F hi_theta_TBT = new H2F("hi_theta_TBT", "hi_theta_TBT", 100, 5.0, 45.0, 20, 5.0, 45.0);  
        hi_theta_TBT.setTitleX("#theta_TBT (deg)"); 
        hi_theta_TBT.setTitleY("#theta_HTCC (deg)");
           
        DataGroup dg_elec = new DataGroup(2,2);
        dg_elec.addDataSet(hi_nphe_all, 0);
        dg_elec.addDataSet(hi_nphe_ele, 0);
        dg_elec.addDataSet(hi_time, 1);
        dg_elec.addDataSet(hi_phi_TBT, 2);
        dg_elec.addDataSet(hi_theta_TBT, 3);
        this.getDataGroup().add(dg_elec,1);
    }
       
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Electrons").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridY(false);
         
        // plotting histos
        this.getAnalysisCanvas().getCanvas("Electrons").cd(0);
        this.getAnalysisCanvas().getCanvas("Electrons").getPad(0).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH1F("hi_nphe_all"));
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH1F("hi_nphe_ele"),"same");
        this.getAnalysisCanvas().getCanvas("Electrons").cd(1);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH1F("hi_time"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(2);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH2F("hi_phi_TBT"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(3);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_TBT"));
        this.getAnalysisCanvas().getCanvas("Electrons").update();
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        DataBank recBankEB   = null;
        DataBank recEvenEB   = null;
        DataBank recDeteEB   = null;
        DataBank recTrajEB   = null;
        DataBank recLTCC     = null;
        DataBank recTBTTrack = null;
        if(event.hasBank("REC::Particle"))  recBankEB = event.getBank("REC::Particle");
        if(event.hasBank("REC::Event"))     recEvenEB = event.getBank("REC::Event");
        if(event.hasBank("REC::Cherenkov")) recDeteEB = event.getBank("REC::Cherenkov");
        if(event.hasBank("REC::Traj"))      recTrajEB = event.getBank("REC::Traj");
        if(event.hasBank("LTCC::clusters"))   recLTCC = event.getBank("LTCC::clusters");
        if(event.hasBank("TimeBasedTrkg::TBTracks")) recTBTTrack = event.getBank("TimeBasedTrkg::TBTracks");
        
        if(recLTCC!=null){
	    int rows = recLTCC.rows();
	    for(int loop = 0; loop < rows; loop++){
                float nphe  = recLTCC.getFloat("nphe",  loop);
                float time  = recLTCC.getFloat("time",  loop);
                float x     = recLTCC.getFloat("x",     loop);
                float y     = recLTCC.getFloat("y",     loop);
                float z     = recLTCC.getFloat("z",     loop);
                float phi   = (recLTCC.getFloat("minPhi",   loop)+recLTCC.getFloat("maxPhi",   loop))/2;
                float theta = (recLTCC.getFloat("minTheta", loop)+recLTCC.getFloat("maxTheta", loop))/2;
                this.getDataGroup().getItem(1).getH1F("hi_nphe_all").fill(nphe);
            }
        }

        if(recBankEB!=null && recEvenEB!=null && recDeteEB!=null && recLTCC!=null && recTrajEB!=null) {
//                    recDeteEB.show(); bank.show();
            double startTime = recEvenEB.getFloat("startTime", 0);
            int nrows = recDeteEB.rows();
	    for(int loop = 0; loop < nrows; loop++){
                int    index = recDeteEB.getShort("index", loop);
                int   pindex = recDeteEB.getShort("pindex", loop);
                int detector = recDeteEB.getByte("detector", loop);
                int pidCode = 0;
                if(recBankEB.getInt("pid", pindex)!=0) pidCode = recBankEB.getInt("pid", pindex);
                else if(recBankEB.getByte("charge", pindex)==-1) pidCode = -211;
                else if(recBankEB.getByte("charge", pindex)==1) pidCode = 211;
                else pidCode = 22;
                Particle recParticle = new Particle(
                                        pidCode,
                                        recBankEB.getFloat("px", pindex),
                                        recBankEB.getFloat("py", pindex),
                                        recBankEB.getFloat("pz", pindex),
                                        recBankEB.getFloat("vx", pindex),
                                        recBankEB.getFloat("vy", pindex),
                                        recBankEB.getFloat("vz", pindex));
                if(recParticle.pid()==11 && detector == DetectorType.LTCC.getDetectorId()) {
                    float nphe  = recLTCC.getFloat("nphe",  index);
                    float time  = recLTCC.getFloat("time",  index);
                    float x     = recLTCC.getFloat("x",     index);
                    float y     = recLTCC.getFloat("y",     index);
                    float z     = recLTCC.getFloat("z",     index);
                    float phi   = (recLTCC.getFloat("minPhi",   index)+recLTCC.getFloat("maxPhi",   index))/2;
                    float theta = (recLTCC.getFloat("minTheta", index)+recLTCC.getFloat("maxTheta", index))/2;
                    double path = 0;
                    for(int j=0; j<recTrajEB.rows(); j++) {
                        if(pindex==recTrajEB.getShort("pindex", j) && recTrajEB.getByte("detector", j)==DetectorType.LTCC.getDetectorId()) {
                            path = recTrajEB.getFloat("path", j);                        
                            this.getDataGroup().getItem(1).getH1F("hi_nphe_ele").fill(nphe);
                            this.getDataGroup().getItem(1).getH1F("hi_time").fill(time-path/29.97-startTime);
                            this.getDataGroup().getItem(1).getH2F("hi_phi_TBT").fill(Math.toDegrees(recParticle.phi()),phi);
                            this.getDataGroup().getItem(1).getH2F("hi_theta_TBT").fill(Math.toDegrees(recParticle.theta()),theta);
                        }
//                            // find corresponding track
//                            double path = 0;
//                            for(int j=0; j<recTBTTrack.rows(); j++) {
//                                if(Math.abs(recParticle.px()-recTBTTrack.getFloat("p0_x", j))<0.1 &&
//                                   Math.abs(recParticle.px()-recTBTTrack.getFloat("p0_x", j))<0.1 &&
//                                   Math.abs(recParticle.px()-recTBTTrack.getFloat("p0_x", j))<0.1) {
//                                    double c3x  = recTBTTrack.getFloat("c3_x",j);
//                                    double c3y  = recTBTTrack.getFloat("c3_y",j);
//                                    double c3z  = recTBTTrack.getFloat("c3_z",j);
//                                    path = recTBTTrack.getFloat("pathlength",j) + Math.sqrt((x-c3x)*(x-c3x)+(y-c3y)*(y-c3y)+(z-c3z)*(z-c3z));
//                    
//                                }
//                            }                    
//                            if(Math.abs(phi-Math.toDegrees(recParticle.phi()))<60 && nFill==0) {
//                                nFill++;
//                            }
    //                        System.out.println(recParticle.phi() + " " + phi);
                        
                    }
                }
            }
    	}       
    }

    @Override
    public void timerUpdate() {

    }


}
