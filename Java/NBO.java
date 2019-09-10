// Copyright (C) 2019 The University of Edinburgh 
// Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
// All Rights Reserved.
// Funded by the Medical Research Council
// https://www.ed.ac.uk/mrc-human-genetics-unit



// naive Bayes object: holding apriori (uniform initialisation) and probability tables (random initialisation)

import java.util.*;
import umontreal.iro.lecuyer.rng.*;


public class NBO {
    public long seed = 12345;
    private MRG32k3a r;
    public int noclasses;
    public int noiter;
    public int nattributes = 4; // in DDD data: 10;
    public static int[] novalues =  {3,3,3,3}; // in DDD data: {2,3, 3,3, 3,3, 3,3, 2,3};
    public double[] apriori;
    public Vector tables;
    public double logl;
    public double[] loglsequence;
    
    public NBO(int nclasses, long seed) {         // create aprori and all tables from a random seed
	setRandom(seed);
	this.noiter = -1;
	this.logl = Double.NEGATIVE_INFINITY;
	this.loglsequence = new double[1];
	this.noclasses = nclasses;
	this.apriori = new double[this.noclasses];
	for(int i=0;  i<this.apriori.length;++i) {
	    this.apriori[i] = 1.0/(double)this.noclasses;
	}
	this.tables = new Vector();
	for(int i=0;  i<nattributes;++i) {
	    double[][] di = new double[this.noclasses][novalues[i]];
	    if(novalues[i]==2) {
		for(int j=0; j<this.noclasses; ++j) {
		    di[j][0] = this.r.nextDouble(); // between 0 and 1
		    di[j][1] = 1-di[j][0];
		}
	    }
	    if(novalues[i]==3) {
		for(int j=0; j<this.noclasses; ++j) {
		    di[j][0] = this.r.nextDouble()/2; // between 0 and 1
		    di[j][1] = this.r.nextDouble()/2; // between 0 and 1
		    di[j][2] = 1-di[j][0]-di[j][1];
		}
	    }
	    this.tables.add(di);
	}
    }

    
    public NBO(int nclasses) {                    // create aprori and all tables with value -1
	this.noiter = -1;
	this.logl = Double.NEGATIVE_INFINITY;
	this.loglsequence = new double[1];
	this.noclasses = nclasses;
	this.apriori = new double[this.noclasses];
	for(int i=0;  i<this.apriori.length;++i) {
	    this.apriori[i] = -1.0;
	}
	this.tables = new Vector();
	for(int i=0;  i<nattributes;++i) {
	    double[][] di = new double[this.noclasses][novalues[i]];
	    if(novalues[i]==2) {
		for(int j=0; j<this.noclasses; ++j) {
		    di[j][0] = -1.0;
		    di[j][1] = -1.0;
		}
	    }
	    if(novalues[i]==3) {
		for(int j=0; j<this.noclasses; ++j) {
		    di[j][0] = -1.0;
		    di[j][1] = -1.0;
		    di[j][2] = -1.0;
		}
	    }
	    this.tables.add(di);
	}
    }
    public void setIter(int a) {
	this.noiter = a;
    }
    public void setLogLikelihood(double d) {
	this.logl = d;
    }
    public void setLogLikelihoodSequence(double[] d) {
	this.loglsequence = new double[d.length];
	for(int i=0;i<this.loglsequence.length;++i) {
	    this.loglsequence[i] = d[i];
	}
    }
    public void setApriori(double [] a) {
	for(int i=0;i<this.noclasses;++i) {
	    this.apriori[i] = a[i];
	}
    }
    public void setTableRow(int i, int j, double[] ijUpdate){            // attribute i class j double[novalues[i]]
	double[][] di = (double[][])this.tables.elementAt(i);
	for(int k=0;k<novalues[i];++k){
	    di[j][k] = ijUpdate[k];
	}
    }
    
    public String toString() {
	String s = new String("Apriori:");
	for(int i=0; i<this.apriori.length;++i) {
	    s=s+this.apriori[i]+"; ";
	}
	String t = new String("Tables:");
	for(int i=0; i<this.tables.size();++i) {
	    t=t+"\n"+(i)+"\n";
	    double[][] di = (double[][])this.tables.elementAt(i);
	    for(int j=0;j<di.length;++j) {
		if(novalues[i]==2) {
		    t=t+di[j][0]+" "+di[j][1];
		}
		if(novalues[i]==3) {
		    t=t+di[j][0]+" "+di[j][1]+" "+di[j][2];
		}
		t=t+"\n";
	    }
	    t=t+"\n---------";
	}
	return(s+"\n"+t+"\niterations:"+this.noiter+"\nlogLikelihood:"+this.logl);
    }
    public String toStringCompact() {
	String s = "";
	for(int i=0; i<this.noclasses;++i) {
	    s += this.apriori[i];                                      // apriori i
	    for(int j=0; j<this.tables.size();++j) {                   // all tables
		double[][] di = (double[][])this.tables.elementAt(j);  // table j
		if(novalues[j]==2) {
		    s += ","+di[i][0]+","+di[i][1];                // row i
		}
		if(novalues[j]==3) {
		    s += ","+di[i][0]+","+di[i][1]+","+di[i][2];
		}
	    }
	    s += "\n";
	}
	return(s);
    }
    public void setRandom(long s) {
	System.out.println("setting random seed:"+s);
	this.seed = s;
	long[] seeds = {s,s,s,s,s,s};
	this.r = null;
	umontreal.iro.lecuyer.rng.MRG32k3a.setPackageSeed(seeds); 
	// see: http://simul.iro.umontreal.ca/ssj/doc/html/umontreal/iro/lecuyer/rng/MRG32k3a.html
	this.r = new MRG32k3a();
    }
    
    public static String toStringArray(int [][] arr) {
	int x = arr.length;
	int y = arr[0].length;
	String result = "";
	for(int i=0;i<x;++i) {
	    for(int j=0;j<y;++j) {
		result += arr[i][j]+" ";
	    }
	    if(i<(x-1)) result += "\n";
	}
	return result;
    }
    public static String toStringArray(double [][] arr) {
	int x = arr.length;
	int y = arr[0].length;
	//if(y>5) { y = 5; }
	String result = "";
	for(int i=0;i<x;++i) {
	    for(int j=0;j<y;++j) {
		result += arr[i][j]+" ";
	    }
	    if(i<(x-1)) result += "\n";
	}
	return result;
    }
    public static String toStringArray(double [] arr) {
	int x = arr.length;
	//if(x>5) { x = 5; }
	String result = "";
	for(int i=0;i<x;++i) {
	    result += arr[i]+"\n";
	}
	return result;
    }
    public static double plusLogLog(double x, double y) {
	if(x==Double.NEGATIVE_INFINITY && y==Double.NEGATIVE_INFINITY)
	    return Double.NEGATIVE_INFINITY;
	
	if(x>y) return (x + Math.log(1 + Math.exp(y-x)));
	else
	    return (y + Math.log(1 + Math.exp(x-y)));
    }
    
    public double loglikeFn(int[][] obs) {
	double []ll = new double[obs.length];
	double llForPZero= Math.log(1e-6);
	for(int i=0;i<obs.length;++i) {                                           // for all obs i
	    double sy = Double.NEGATIVE_INFINITY;                                // sum for all classes sy
	    for(int y=0;y<this.noclasses;++y) {                                   // for all classes y
		double piy = Math.log(this.apriori[y]);                           // piy is apriori for y
		for(int j=0;j<this.tables.size();++j) {                           // for all attributes j
		    double[][] atj = (double[][])this.tables.elementAt(j);        // the jth table atj
		    if(obs[i][j]>=0) {                                            // if obs value not NA (-1) :index into table column
			//System.out.println("obs:"+i+"; class:"+y+"; attr:"+j+"; obs val:"+obs[i][j]+"; p val:"+atj[y][obs[i][j]]+"\n");
			piy += Math.max(llForPZero,Math.log(atj[y][obs[i][j]]));  // select row y (table j) column obs[i][j] 
			//System.out.println("piy"+piy+"\n");
		    } else {
			//System.out.println("NA:: obs:"+i+"; class:"+y+"; attr:"+j+"; obs val:"+obs[i][j]+"; p val: 1/"+this.novalues[j]+"\n");
			piy += Math.log(1.0/(double)this.novalues[j]);                      // for NA use log(1/no.values)
			//System.out.println("piy"+piy+"\n");
		    }
		}
		sy = plusLogLog(sy,piy);                                        // add (logged) likelihood values for all attributes
	    }
	    ll[i] = sy;                                                          // save sum for obs  
	}
	double sm = 0.0;
	for(int i=0;i<obs.length;++i) {
	    if(ll[i]!=-Double.NEGATIVE_INFINITY) { sm += ll[i]; }  // ll is sum log l (when not -Inf)
	    //if(i<10) { System.out.println("ll: "+i+":"+ll[i]);  }            
	}
	return(sm);
    }

    

    public static double productOverAttributes(NBO bo, int j, int [] obsi) {
	double r = 1;
	for(int i=0;i<bo.nattributes;++i) {
	    if(obsi[i]!=-1) {
		double[][] di = (double[][])bo.tables.elementAt(i);
		r = r * di[j][obsi[i]];
	    }
	}
	return(r*bo.apriori[j]);
    }
    
    public static NBO em_nBFn(int[][] obs, NBO bo, int maxiter, double minDeltaLL) {
	double lli = 0;
	double ll = bo.loglikeFn(obs);                                 // ll of the data given initial (random) tables
	int noObs = obs.length;
	int noAttr = obs[0].length;
	int noClasses = bo.noclasses;
	int itn = 0;
	double[] logls =new double[maxiter];
	while(itn<maxiter) {
	    if(itn>0) { ll = lli; }
	    double[][] delta = new double[noClasses][noObs];            // responsibility array classes * obs
	    for(int i=0;i<noObs;++i) {
		for(int j=0;j<noClasses;++j) {
		    delta[j][i] = productOverAttributes(bo,j,obs[i]);   // for class j, for all attributes, get product pr value=obs value * apriori[j]
		}
	    }
	    //System.out.println("delta:"+toStringArray(delta));
	                                                                // set responsibilities for obs (columns) to sum to 1
	    double[] sumYs = new double[noObs];                         // sum responsibilities per obs
	    for(int i=0;i<noObs;++i) {
		double si = 0;
		for(int j=0;j<noClasses;++j) {
		    si += delta[j][i];
		}
		sumYs[i] = si;
	    }
	    //System.out.println("sumYs:"+toStringArray(sumYs));
	    for(int i=0;i<noObs;++i) {                                   
		if(sumYs[i]!=0) {                                        // divide column i for all j [j][i] by sumYs[i]
		    for(int j=0;j<noClasses;++j) {
			delta[j][i] = delta[j][i]/sumYs[i];
		    }
		} else {
		    System.out.println("Warning: sumY=0; ");
		    for(int j=0;j<noClasses;++j) {
			delta[j][i] = 1.0/(double)noClasses;
		    }
		}
	    }
	    //System.out.println("delta:"+toStringArray(delta));
	    double[] aprioriUpdate = new double[noClasses];            // update apriori
	    for(int j=0;j<noClasses;++j) {
		double sj = 0;
		for(int i=0;i<noObs;++i) {
		    sj += delta[j][i];
		}
		aprioriUpdate[j] = sj/(double)noObs;
	    }
	    double sumAprioriUpdate = 0;
	    for(int j=0;j<noClasses;++j) {
		sumAprioriUpdate += aprioriUpdate[j];
	    }
	    if(sumAprioriUpdate!=1.0) {
		for(int j=0;j<noClasses;++j) {
		    aprioriUpdate[j] = aprioriUpdate[j]/sumAprioriUpdate;
		}
	    }
	    //System.out.println("aprioriUpdate:"+toStringArray(aprioriUpdate));
	    NBO boUpdate = new NBO(noClasses);                     // new NB object
	    boUpdate.setApriori(aprioriUpdate);                    // set apriori
	    for(int i=0;i<noAttr;++i) {                            // for all attributes i
		int[] valsi = new int [novalues[i]];                           // valsi the allowable values in obs[][i]
		if(novalues[i]==2) {
		    valsi[0] = 0; valsi[1] = 1; 
		} else {
		    valsi[0] = 0; valsi[1] = 1; valsi[2] = 2; 
		}
		for(int j=0;j<noClasses;++j) {                     // for all classes j
		    double[] ijUpdate = new double[novalues[i]];   // update for attr i class j
		    for(int k=0;k<valsi.length;++k) {              // for all k in valsi
			double sumjk = 0;                          // sum delta[j][] for obs matching valsi[k]
			double sumjkall = 0;                       // sum all delta[j][] 
			for(int l=0;l<noObs;++l) {
			    if(obs[l][i]==valsi[k]) {
				sumjk += delta[j][l];
			    }
			    sumjkall += delta[j][l];
			}
			ijUpdate[k] = sumjk/sumjkall;              // update is ratio of sum matching to sum all 		    
		    }
		    double smIJUpdate = 0;
		    for(int k=0;k<valsi.length;++k) { smIJUpdate += ijUpdate[k]; }
		    if(smIJUpdate==0) {
			for(int k=0;k<valsi.length;++k) { ijUpdate[k] = 1.0/(double)novalues[i]; }
		    } else if(smIJUpdate!=0){
			for(int k=0;k<valsi.length;++k) { ijUpdate[k] = ijUpdate[k]/smIJUpdate; }
		    }
		    boUpdate.setTableRow(i,j,ijUpdate);           // attribute i class j double[novalues[i]]
		}
	    }
	    //System.out.println(boUpdate.toString());
	    bo = boUpdate;
	    lli = bo.loglikeFn(obs);
	    bo.setLogLikelihood(lli);
	    logls[itn] = lli;
	    bo.setLogLikelihoodSequence(logls);
	    //System.out.println("lli:"+lli);
	    if(itn>50 && Math.abs(lli-ll)<minDeltaLL) { bo.setIter(itn); return(bo); }
	    itn += 1;
	}
	if(Math.abs(lli-ll)>=minDeltaLL) { System.out.println("Warning: min delta LL exceeded at max iter "+(lli-ll)); }
	bo.setIter(itn);
	return(bo);
    }

     public static void main(String[] args) {
	 int noclasses;
	 int seed;
	 String datafile;
	 String runid;
	 if(args.length==3) {
	     seed = Integer.parseInt(args[0]);
	     noclasses = Integer.parseInt(args[1]);
	     datafile = new String(args[2]);
	     runid = datafile.substring(0,datafile.indexOf("."));
	     System.out.println("running EM NBayes with seed="+seed+" k="+noclasses);
	     Read ro = new Read(datafile);
	     Vector v  = ro.readData(ro.datafile);
	     int[][] obs = Read.toArray(v);
	     System.out.println("read data:"+obs.length+"*"+obs[0].length+"\n");
	     NBO r = new NBO(noclasses,seed);
	     long t1 = System.currentTimeMillis();
	     NBO s = em_nBFn(obs, r, 1000, 1e-3);  // 1000 for 30 classes
	     long t2 = System.currentTimeMillis();
	     Read.setFileContent(".","EMNB_"+seed+"_"+noclasses+"_"+runid+"_logls.txt",toStringArray(s.loglsequence));
	     Read.setFileContent(".","EMNB_"+seed+"_"+noclasses+"_"+runid+"_model.txt",s.toStringCompact());
	     double ll = s.logl;
	     System.out.println("time [0]: "+(t2-t1)+"(ms)"+" ll:"+ll);
	     for(int i=1;i<1000;++i) {                
		 r = new NBO(noclasses,seed+(3*i)); 
		 t1 = System.currentTimeMillis();
		 s = em_nBFn(obs, r, 1000, 1e-3);
		 t2 = System.currentTimeMillis();
		 if(s.logl>ll) {
		     Read.setFileContent(".","EMNB_"+seed+"_"+noclasses+"_"+runid+"_logls.txt",toStringArray(s.loglsequence));
		     Read.setFileContent(".","EMNB_"+seed+"_"+noclasses+"_"+runid+"_model.txt",s.toStringCompact());
		     ll = s.logl;
		     System.out.println("time ["+i+"]: "+(t2-t1)+"(ms)"+" ll:"+s.logl+" **best**");		     
		 } else {		     
		     System.out.println("time ["+i+"]: "+(t2-t1)+"(ms)"+" ll:"+s.logl);
		 }
	     }
	 }
     }
		

}
