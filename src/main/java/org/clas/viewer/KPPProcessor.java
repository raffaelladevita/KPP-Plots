package org.clas.viewer;

import java.util.Arrays;
import org.clas.analysis.HBTmonitor;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.Collectors;
import org.clas.analysis.CNDmonitor;
import org.clas.analysis.CTOFmonitor;
import org.clas.analysis.CVTmonitor;
import org.clas.analysis.EBHBmonitor;
import org.clas.analysis.EBmonitor;
import org.clas.analysis.ECmonitor;
import org.clas.analysis.ELASTICmonitor;
import org.clas.analysis.ELmonitor;
import org.clas.analysis.EPI0monitor;
import org.clas.analysis.EPIPLUSmonitor;
import org.clas.analysis.FMTmonitor;
import org.clas.analysis.FTOFmonitor;
import org.clas.analysis.FTmonitor;
import org.clas.analysis.HTCCmonitor;
import org.clas.analysis.KINEmonitor;
import org.clas.analysis.LTCCmonitor;
import org.clas.analysis.TBTmonitor;
import org.clas.analysis.TIMEmonitor;

import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

/**
 *
 * @author ziegler
 */
public class KPPProcessor {

    String prefix;
    String flist;
    CodaEventDecoder decoder = new CodaEventDecoder();
    DetectorEventDecoder detectorDecoder = new DetectorEventDecoder();

    TreeMap<String, List<H2F>> histos = new TreeMap<String, List<H2F>>();

    private int runNumber = 0;

    // detector monitors
    AnalysisMonitor[] monitors = {
        new HBTmonitor("HBT"),
        new TBTmonitor("TBT"),
        new CVTmonitor("CVT"),
        new ELmonitor("ELECTRONS"),
        new ECmonitor("EC"),
        new FTOFmonitor("FTOF"),
        new HTCCmonitor("HTCC"),
        new LTCCmonitor("LTCC"),
        new CTOFmonitor("CTOF"),
        new CNDmonitor("CND"),
        new FTmonitor("FT"),
        new FMTmonitor("FMT"),
        new EBHBmonitor("EBHB"),
        new EBmonitor("EB"),
        new KINEmonitor("KINEMATICS"),
        new ELASTICmonitor("ELASTIC"),
        new EPIPLUSmonitor("EPIPLUS"),
        new EPI0monitor("EPI0"),
        new TIMEmonitor("TIME")
    };

    public KPPProcessor(String prefix, String list) {
        this.prefix = prefix;
        this.flist = list;
    }

    private int getRunNumber(DataEvent event) {
        int rNum = this.runNumber;
        DataBank bank = event.getBank("RUN::config");
        if (bank != null) {
            rNum = bank.getInt("run", 0);
        }
        return rNum;
    }

    public boolean processDataEvent(DataEvent event) {
        if (event != null) {
            if (this.runNumber != this.getRunNumber(event)) {
                if (this.runNumber != 0) {
                    saveHistosToFile();
                    resetEventListener();
                }
                this.runNumber = this.getRunNumber(event);
            }

            for (int k = 0; k < this.monitors.length; k++) {
                this.monitors[k].dataEventAction(event);
            }
        }
        return true;
    }

    public void resetEventListener() {
        for (int k = 0; k < this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
        }
    }

    public void saveHistosToFile() {
        TDirectory dir = new TDirectory();
        for (int k = 0; k < this.monitors.length; k++) {
            this.monitors[k].writeDataGroup(dir);
        }
        dir.writeFile(prefix+"_"+flist+".hipo");
    }

    public boolean init() {
        return true;
    }

    public static void main(String[] args) throws Exception {
        String flist = Arrays.stream(args)
                .map(s -> s.substring(s.lastIndexOf("/") + 1))
                .collect(Collectors.joining(";"));
        
        KPPProcessor kppProc = new KPPProcessor("kpp_histos", flist);

        for (int iarg = 0; iarg < args.length; iarg++) {
            HipoDataSource reader = new HipoDataSource();
            reader.open(args[iarg]);
            while (reader.hasEvent() == true) {
                DataEvent event = reader.getNextEvent();
                kppProc.processDataEvent(event);
            }
        }
        kppProc.saveHistosToFile();
    }
}
