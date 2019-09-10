// Copyright (C) 2019 The University of Edinburgh 
// Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
// All Rights Reserved.
// Funded by the Medical Research Council
// https://www.ed.ac.uk/mrc-human-genetics-unit

import java.io.*;
import java.util.*;


public class Read {
    String datafile;
    public Read(String f) {
	this.datafile = f;
    }
    public Read() {
	this.datafile = "initModelsData.csv";
    }
    
    
    public static Vector readData(String filename) {
	String r;
	Vector v = new Vector();
	DO d;
	try {
	    BufferedReader in = new BufferedReader(new FileReader(filename));
	    while((r=in.readLine())!=null) {
		//System.out.println(r);
		v.addElement(d = new DO(r,NBO.novalues));
		//System.out.println(d.toString());
	    }
	} catch(Exception e) {e.printStackTrace(); System.out.println("Read.java - exception reading file");}
	return(v);
    }

    public static int[][] toArray(Vector v)  {
	int[][] a = new int[v.size()][NBO.novalues.length];
	for(int i=0;i<v.size();++i) {
	    a[i] = ((DO)v.elementAt(i)).d;
	}
	return(a);
    }

    public static void setFileContent(String dir, String filename, String content) {
	setFileContent(dir,filename,content,true);
    }
    public static void setFileContent(String dir, String filename, String content,boolean report) {
	if((filename == null) || (filename.length() == 0)) return;
	File f;
	FileWriter out = null;
	try {
	    f = new File(dir,filename);
	    out = new FileWriter(f);
	    // very simplistic approach !
	    int size = content.length();
	    out.write(content,0,size);
	}
	catch (IOException e) {
	    System.out.println("file write error: "+dir+" "+filename);
	    if(report) System.out.println(content);
	    return; }
	finally {
	 try {
	     if(out != null) out.close();} catch (IOException e) {}
	}
    }
    
    public static void main(String[] args) {
	Read r = new Read();
	Vector v  = r.readData(r.datafile);
    }
    
}
