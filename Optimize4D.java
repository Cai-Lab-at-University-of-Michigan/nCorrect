import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URISyntaxException;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.CompositeImage;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.ChannelSplitter;

/**
 * @author Sara Azzouz <sazzouz@umich.edu>
 */

public class Optimize4D {
	
	private static int maxX;
	private static int maxY;
	private static int maxZ;
	private static int nChannel;
	private static ImageStack[] istk;
	private static int typeId;
	private static double[] colorRef;
	private static ArrayList<int[]> oPathN = new ArrayList<>();
	private static ArrayList<double[]> pPath = new ArrayList<>();
	
	private static double aweight = 1.00;
	private static double bweight = 0.85;
	private static double cweight = 3.75;
	private static double cossimT = 0.90;
	private static double backgrT = 0.47;
	private static double minIntT = 10000;
	private static double maxIntT = 85000;
	private static double minPerC = 0.05;
	private static double maxPerC = 0.30;
	private static int maxR = 12;
	
	public static void main(String[] args) throws IOException, URISyntaxException{
		if(args.length!=13)
			System.out.println("Error: Invalid number of parameters!");
 		else {
 			String allAddress = args[0].replaceAll("\\\\","/")+"/";
			String imgAddress = args[1].replaceAll("\\\\","/")+"/";
			String outAddress = args[2].replaceAll("\\\\","/")+"/";
			aweight = Double.parseDouble(args[3]);
			bweight = Double.parseDouble(args[4]);
			cweight = Double.parseDouble(args[5]);
			cossimT = Double.parseDouble(args[6]);
			backgrT = Double.parseDouble(args[7]);
			minIntT = Double.parseDouble(args[8]);
			maxIntT = Double.parseDouble(args[9]);
			minPerC = Double.parseDouble(args[10]);
			maxPerC = Double.parseDouble(args[11]);
			maxR    = Integer.parseInt(args[12]);
			
			MyFileInfoVirtualStack example = new MyFileInfoVirtualStack(imgAddress,false);
			ImagePlus imp  = example.getImage();
			maxX = imp.getHeight();
			maxY = imp.getWidth();
			maxZ = imp.getNSlices();
			nChannel = imp.getNChannels();
			int maxValue = (int)(Math.pow(2, imp.getBitDepth())-1);
			
			istk = new ImageStack[nChannel];
			for(int i = 1; i <= nChannel; ++i) {
				istk[i-1]= ChannelSplitter.getChannel(new ImagePlus(imgAddress),i);
			}
			
			ArrayList<String> addresses = readAddressFile(allAddress);
			PrintStream orig = System.out; // for debugging only
			PrintStream stats = new PrintStream(new FileOutputStream(outAddress+"stats.csv"));
			System.setOut(stats);
			System.out.println("SWC File Name"+", "+"Original Path Length"+", "+"Corrected Path Length"+", "+"Pathfinding Time"+", "+"Total Run Time");
			
			for(int swcInd = 0; swcInd < addresses.size(); ++swcInd) {
				long startNeurTime = System.currentTimeMillis();
				String swcAddress = addresses.get(swcInd);
				ArrayList<double[]> oPath = readSWCFile(swcAddress);
				ArrayList<Integer> pointPairs = findSEPointPairs(oPath);
				
				ImageStack pStack = ImageStack.create(maxY, maxX, maxZ, imp.getBitDepth());
				ImageStack oStack = ImageStack.create(maxY, maxX, maxZ, imp.getBitDepth());
				ImageStack raInfo = ImageStack.create(maxY, maxX, maxZ, imp.getBitDepth());

				System.setOut(orig); // for debugging only
				long startPathTime = System.currentTimeMillis();
				while(!pointPairs.isEmpty()){	
					int startInd = pointPairs.remove(0);
					int endInd   = pointPairs.remove(0);
					
					double[] pt = oPath.get(startInd);			
					if(pt[6] != -1) {
						int[] endPt = new int[]{(int)pt[3], (int)pt[2], (int)pt[4]};	
						int minIndex = 0;
						double minDist = Double.POSITIVE_INFINITY;
						for(int p = 0; p < pPath.size(); ++p) {
							pt = pPath.get(p);
							double tempDist = (pt[3]-endPt[0])*(pt[3]-endPt[0])+(pt[2]-endPt[1])*(pt[2]-endPt[1])+(pt[4]-1-endPt[2])*(pt[4]-1-endPt[2]);
							if(tempDist <= minDist) {
								minIndex = p;
								minDist = tempDist;
							}
						}
						pt = pPath.get(minIndex);
						int[] strPt = new int[]{(int) pt[3], (int) pt[2], (int) pt[4]-1};
						
						Pathfinder4D pf = new Pathfinder4D(new Grid4D(maxX, maxY, maxZ));
						ArrayList<int[]> path = pf.findPath(strPt,endPt);
						
						for(int i = 0; i <= path.size()-1; ++i) {
							int[] p = path.get(i);
							oPathN.add(new int[] {p[0], p[1], p[2], 0});
						}
					}
					
					for(int i = startInd; i <= endInd; ++i) {
						pt = oPath.get(i);
						oPathN.add(new int[] {(int)pt[3], (int)pt[2], (int)pt[4], 0});
					}
					colorRef = avgColor(0, oPathN.size(), false);
					
					double totPenalty = 0;
					for(int i = 0; i < oPathN.size(); ++i) {
						double[] bestRad = getBestRadius(oPathN.get(i), i);
						oPathN.get(i)[3] = (int)bestRad[0];
						totPenalty += bestRad[1];
					}
					System.out.println(swcInd+" START: "+totPenalty);
					
					double pathPenalty = totPenalty;
					int numConst = 0;
					while(numConst == 0) {
						if(oPathN.size()<=2)
							break;
						for(int i=1; i < oPathN.size()-1; ++i) {
							ArrayList<int[]> sharedNeighbors = getSharedNeighbors(oPathN.get(i-1),oPathN.get(i+1),oPathN.get(i));
							int ind = -1;
							double[] pair = new double[] {0, Double.POSITIVE_INFINITY};
							for(int j = 0; j < sharedNeighbors.size(); ++j) {
								double[] tempPair = getBestRadius(sharedNeighbors.get(j), i);			
								if(tempPair[1] < pair[1]) {
									pair = tempPair;
									ind = j;
								}
							}
							if(ind != -1) {
								int[] newPt = sharedNeighbors.get(ind);
								newPt[3] = (int) pair[0];
								double oldPenalty = calculateColorBonus(oPathN.get(i), i);
								pathPenalty -= oldPenalty-pair[1];
								oPathN.set(i, newPt);
							}
						}		
						System.out.println(swcInd+"      : "+pathPenalty);
						if(pathPenalty < totPenalty) {
							totPenalty = pathPenalty;
							numConst = 0;
						}
						else
							numConst++;
					}
					System.out.println(swcInd+" AFTER: "+pathPenalty);
					appendNewPath(oPathN);
					oPathN.clear();
				}
				long stopPathTime = System.currentTimeMillis();
				long elapPathTime = stopPathTime-startPathTime;
				
				String swcFileName = swcAddress.substring(swcAddress.lastIndexOf("/")+1,swcAddress.lastIndexOf("."));	
				PrintStream out = new PrintStream(new FileOutputStream(outAddress+"corrected_path_"+swcFileName+".swc"));
				System.setOut(out);
				
				for(int i = 0; i < pPath.size(); ++i) {
					double[] point = pPath.get(i);
					System.out.println((int)point[0]+" "+(int)point[1]+" "+point[2]+" "+point[3]+" "+point[4]+" "+point[5]+" "+(int)point[6]);
					pStack.getProcessor((int)point[4]).set((int)point[2], (int)(point[3]), maxValue); 
				}
				
				for(int i = 0; i < oPath.size(); ++i) {
					double[] point = oPath.get(i);
					oStack.getProcessor((int)point[4]+1).set((int)point[2], (int)point[3], maxValue);
				}
				
				for(int i = 0; i < pPath.size(); ++i) {
					double[] point = pPath.get(i);
					int[] centerPoint = new int[] {(int)point[3], (int)point[2], (int)point[4]-1, (int)point[5]};
					int r = centerPoint[3];
					for(int z = centerPoint[2]-r; z <= centerPoint[2]+r; ++z) {
			    		for(int x = centerPoint[0]-r; x <= centerPoint[0]+r; ++x) {
			    			for(int y = centerPoint[1]-r; y <= centerPoint[1]+r; ++y) {
								int dx = x-centerPoint[0];
								int dy = y-centerPoint[1];
								int dz = z-centerPoint[2];
								int distSquare = dx*dx+dy*dy+dz*dz;
			    				if(isValid(x,y,z) && distSquare <= r*r)
			    					raInfo.getProcessor(z+1).set(y, x, maxValue);
			    			}
			    		}
			    	}			
				}
	
				ImagePlus[] net = new ImagePlus[nChannel+3];
				for(int i = 0; i < nChannel+3; ++i) {
					net[i] = new ImagePlus();
					if(i < nChannel)
						net[i].setStack(istk[i]);
				}
				net[nChannel].setStack(pStack);
				net[nChannel+1].setStack(oStack);
				net[nChannel+2].setStack(raInfo);
	
				CompositeImage c = new CompositeImage(RGBStackMerge.mergeChannels(net, true));
				c.setDimensions(nChannel+3, maxZ, 1);
				c.setDisplayMode(IJ.COMPOSITE);
				IJ.save(c, outAddress+"corrected_imag_"+swcFileName+".tif");
	
				ZProjector zProjector = new ZProjector(net[nChannel]);
				zProjector.setMethod(ZProjector.MAX_METHOD);
				zProjector.doProjection(true);
				IJ.save(zProjector.getProjection(), outAddress+"corrected_proj_"+swcFileName+".tif");
								
				long stopNeurTime = System.currentTimeMillis();
				long elapNeurTime = stopNeurTime-startNeurTime;

				System.setOut(stats);
				System.out.println(swcFileName+", "+oPath.size()+", "+pPath.size()+", "+elapPathTime+", "+elapNeurTime);
				pPath.clear();
			}
 		}
	}
	
	private static ArrayList<String> readAddressFile(String address) throws FileNotFoundException{
		ArrayList<String> addresses = new ArrayList<>();
		Scanner s = new Scanner(new File(address));
		while(s.hasNextLine()) {
			addresses.add(s.nextLine());
		}
		s.close();		
		return addresses;
	}	
	private static ArrayList<double[]> readSWCFile(String address) throws FileNotFoundException{
		ArrayList<double[]> original = new ArrayList<>();
		Scanner s = new Scanner(new File(address));
		while(s.hasNextLine()) {
			String input = s.nextLine();
			Scanner r = new Scanner(input);
			if(!input.isEmpty() && !input.substring(0,1).contentEquals("#"))
				original.add(new double[] {r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble()-1, r.nextDouble(), r.nextDouble()});
			r.close();
		}
		s.close();
		typeId = (int) original.get(0)[1];
		return original;
	}
	private static void appendNewPath(ArrayList<int[]> path) {
		for(int i = 0; i < path.size(); ++i) {
			int[] point = path.get(i);
			if(i == 0) {
				int pPathInd = -1;
				for(int j = 0; j < pPath.size(); ++j) {
					int[] pPoint = new int[] {(int) pPath.get(j)[3], (int) pPath.get(j)[2], (int) pPath.get(j)[4]};
					if(pPoint[0] == point[0] && pPoint[1] == point[1] && pPoint[2] == point[2]+1) {
						pPathInd = j;
						break;
					}
				}
				if(pPathInd == -1)					
					pPath.add(new double[] {1+pPath.size(), typeId, point[1], point[0], point[2]+1, point[3], -1});
				else {
					double[] refPoint = pPath.get(pPathInd);
					refPoint[1] = 5;
					pPath.set(pPathInd, refPoint);
					point = path.get(++i);
					pPath.add(new double[] {1+pPath.size(), typeId, point[1], point[0], point[2]+1, point[3], pPathInd+1});
				}
			}
			else if(i != path.size()-1)
				pPath.add(new double[] {1+pPath.size(), typeId, point[1], point[0], point[2]+1, point[3], pPath.size()});
			else
				pPath.add(new double[] {1+pPath.size(), 6, point[1], point[0], point[2]+1, point[3], pPath.size()});
		}
	}	
	private static ArrayList<Integer> findSEPointPairs(ArrayList<double[]> list){
		ArrayList<Integer> pointPairs = new ArrayList<>();
		int i = 0;
		while(i < list.size()) {
			double[] pt = list.get(i);
			pointPairs.add(i);
			while(i < list.size() && pt[1] != 6)
				pt = list.get(++i);
			pointPairs.add(i++);
		}
		return pointPairs;
	}
	private static boolean isValid(int xpos, int ypos, int zpos){
        return xpos >= 0 && ypos >= 0 && zpos >= 0 && xpos < maxX && ypos < maxY && zpos < maxZ;
    }
	private static double[] avgColor(int startInd, int endInd, boolean isAng) {
		double[] colorRef = new double[nChannel];
		double[] colorPnt = new double[nChannel];
		for(int i = startInd; i < endInd; ++i) {
			int[] pt = {oPathN.get(i)[0], oPathN.get(i)[1], oPathN.get(i)[2]};			
			for(int channel = 0; channel < nChannel; ++channel)
				colorPnt[channel] = istk[channel].getVoxel(pt[1], pt[0], pt[2]);
			if(isAng)
				colorPnt = norm(colorPnt);
			for(int channel = 0; channel < nChannel; ++channel) {
				colorRef[channel] += colorPnt[channel]/(endInd-startInd+1.0);
			}
		}
		return colorRef;
	}
	private static int sum(double[] color) {
		int sum = 0;
		for(int i = 0; i < nChannel; ++i)
			sum += color[i];
		return sum;
	}
	private static double[] norm(double[] vec) {
		double mag = 0;
		for(int i = 0; i < nChannel; ++i)
			mag += vec[i]*vec[i];
		mag = Math.sqrt(mag);
		for(int i = 0; i < nChannel; ++i) {
			vec[i] *= 1.0/mag;
		}
		return vec;
	}
	public static double[] getBestRadius(int[] pt, int nearestInd) {
		int rpos = pt[3];
		double minPenalty = Double.POSITIVE_INFINITY;
		for (int r = 1; r < maxR; ++r) {
			pt[3] = r;
			double tempPenalty = calculateColorBonus(pt, nearestInd);
			if(tempPenalty <= minPenalty) {
				minPenalty = tempPenalty;
				rpos = r;
			}
		}
		return new double[] {rpos, minPenalty};
	}
	public static double calculateColorBonus(int[] pos, int nearestInd) {
		int r = pos[pos.length-1];
		ArrayList<int[]> indexes = new ArrayList<int[]>();
		for(int z = pos[2]-r; z <= pos[2]+r; ++z) {
			for(int x = pos[0]-r; x <= pos[0]+r; ++x) {
				for(int y = pos[1]-r; y <= pos[1]+r; ++y) {
					int dx = x-pos[0];
					int dy = y-pos[1];
					int dz = z-pos[2];
					int distSquare = dx*dx+dy*dy+dz*dz;
					if(isValid(x,y,z) && distSquare < r*r) {
						indexes.add(new int[] {x,y,z});
					}						
				}
			}
		}
		
		int refColSum = sum(colorRef);
		
		double percentC = 0;
		double m = (maxPerC-minPerC)/(minIntT-maxIntT);
		double b = maxPerC-m*minIntT;
		if(refColSum <= minIntT)
			percentC = maxPerC;
		else if(refColSum > minIntT && refColSum < maxIntT)
			percentC = m*refColSum+b;
		else
			percentC = minPerC;
		
		double percentIntB = 0;		
		double percentColB = 0;
		for(int i = 0; i < indexes.size(); ++i) {
			int[] index = indexes.get(i);
			double[] nodeColor = new double[nChannel];
			for(int channel = 0; channel < nChannel; ++channel)
				nodeColor[channel] = istk[channel].getVoxel(index[1],index[0],index[2]);	
			double dotProd = 0; double magNColor = 0; double magRColor = 0;
			for(int l = 0; l < nChannel; ++l) {
				dotProd   += nodeColor[l]*colorRef[l];
				magNColor += nodeColor[l]*nodeColor[l];
				magRColor += colorRef[l]*colorRef[l];
			}
			boolean angCond = (dotProd/Math.sqrt(magNColor*magRColor) <= cossimT) ? false : true;
			boolean intCond = (sum(nodeColor) <= percentC*refColSum || 
							   sum(nodeColor) >= 3.0*refColSum) ? false : true;	
			percentIntB += intCond ? 0 : 1;
			percentColB += (angCond || intCond) ? 0 : 1;
		} 
		percentIntB *= 1.0/(indexes.size()*1.0);
		percentColB *= 1.0/(indexes.size()*1.0);   

		if (percentIntB >= backgrT && r>1)
			percentIntB = 100.0;
		double penalty = aweight*percentIntB+bweight*percentColB+cweight/(r*r);
		return penalty;
	}
	public static ArrayList<int[]> getSharedNeighbors(int[] p1, int[] p2, int[] priority){
		ArrayList<int[]> neighbors = new ArrayList<>();
		HashMap<String, int[]> p1Neighbors = new HashMap<>();		
		for(int x = -1; x <= 1; ++x){
			for(int y = -1; y <= 1; ++y){
				for(int z = -1; z <= 1; ++z){
					if(x == 0 && y == 0 && z == 0) {continue;}
					int xpos = p1[0]+x;
					int ypos = p1[1]+y;
					int zpos = p1[2]+z;
					if(isValid(xpos,ypos,zpos))
						p1Neighbors.put(xpos+" "+ypos+" "+zpos, new int[]{xpos,ypos,zpos,0});	
				}
			}
		}
		for(int x = -1; x <= 1; ++x){
			for(int y = -1; y <= 1; ++y){
				for(int z = -1; z <= 1; ++z){
					if(x == 0 && y == 0 && z == 0) {continue;}
					int xpos = p2[0]+x;
					int ypos = p2[1]+y;
					int zpos = p2[2]+z;
					if(xpos == priority[0] && ypos == priority[1] && zpos == priority[2])
						neighbors.add(0, new int[] {xpos, ypos, zpos, 0});
					else if(isValid(xpos,ypos,zpos) && p1Neighbors.containsKey(xpos+" "+ypos+" "+zpos))
						neighbors.add(new int[] {xpos, ypos, zpos, 0});	
				}
			}
		}
		return neighbors;
	}    
}