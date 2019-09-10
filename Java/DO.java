// Copyright (C) 2019 The University of Edinburgh 
// Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
// All Rights Reserved.
// Funded by the Medical Research Council
// https://www.ed.ac.uk/mrc-human-genetics-unit

import java.util.*;

public class DO {
    public int[] d;
    public DO(String s,int[] novalues) {
	this.d = new int[novalues.length];
	StringTokenizer st = new StringTokenizer(s,",");
	int i=0;
	String t;
	while (st.hasMoreTokens()) {
	    t = st.nextToken();
	    if(t.equals("NA")) { this.d[i] = -1;
	    } else if(t.equals("-1")) {  this.d[i] = 0;
	    } else if(t.equals("0")) {  this.d[i] = 1; if(novalues[i]==2) {  System.out.println("**[1] error unrecognised input**"+t); }
	    } else if(t.equals("1")) {  if(novalues[i]==3) { this.d[i] = 2; } else if(novalues[i]==2) { this.d[i] = 1; } else { System.out.println("**[2] error unrecognised input**"+t); }
	    } else { System.out.println("**[3] error unrecognised input**"+t); }
	    i = i+1;
	}
    }
    
    public String toString() {
	String s = new String();
	for(int i=0; i<d.length;++i) {
	    s=s+d[i];
	}
	return(s);
    }

}
