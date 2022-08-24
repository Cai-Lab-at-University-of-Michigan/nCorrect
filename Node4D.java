/**
 * @author Sara Azzouz <sazzouz@umich.edu>
 */

public class Node4D implements Comparable<Node4D>{
    
    private int[]  position;
    private Node4D parent;       
    private double gCost;
    private double hCost;   
  
    public Node4D(int[] position){
        this.position = position;
    } 
   
    public int[] getPosition(){
        return position;
    }
    public Node4D getParent(){
        return parent;
    }
    public double getGCost(){
        return gCost;
    }
    public double getHCost(){
        return hCost;
    }
    public double getFCost(){
        return gCost + hCost;
    }

    public void setParent(Node4D parent){
        this.parent = parent;
    }
    public void setGCost(double gCost){
        this.gCost = gCost;
    }
    public void setHCost(double hCost){
        this.hCost = hCost;
    }

    @Override
    public int compareTo(Node4D n) {
        int compare = ((Double)getFCost()).compareTo(n.getFCost());
        if(compare == 0){
            compare = ((Double)getHCost()).compareTo(n.getHCost());
        }
        return -compare;
    }     
}