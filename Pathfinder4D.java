import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.PriorityQueue;

/**
 * @author Sara Azzouz <sazzouz@umich.edu>
 */

public class Pathfinder4D {
    
    private Grid4D gr;

    public Pathfinder4D(Grid4D gr){
        this.gr = gr;
    }

    public ArrayList<int[]> findPath(int[] start, int[] end){
    	if(gr.isValid(start[0], start[1], start[2]) && gr.isValid(end[0], end[1], end[2])){
            PriorityQueue<Node4D> openSet = new PriorityQueue<>(Collections.reverseOrder());
            HashSet<Node4D> closedSet = new HashSet<>();
            
            Node4D endNode = gr.getNode(end[0], end[1], end[2]);
            Node4D startNode = gr.getNode(start[0], start[1], start[2]);
            
            openSet.add(gr.getNode(start[0], start[1], start[2]));
            
            while(!openSet.isEmpty()){
                Node4D currentNode = openSet.poll();
                closedSet.add(currentNode);
                
                if(currentNode == endNode){
                    return retrace(startNode,endNode);                    
                }
                for(Node4D neighbor : gr.getNeighbors(currentNode)){
                    if(closedSet.contains(neighbor)) {continue;}
                                        
                    double newMoveCost2Neighbor = currentNode.getGCost()+gr.getDistance(currentNode, neighbor);
                    if(newMoveCost2Neighbor < neighbor.getGCost() || !openSet.contains(neighbor)) {
                        neighbor.setGCost(newMoveCost2Neighbor);
                        neighbor.setHCost(gr.getDistance(neighbor, endNode));
                        neighbor.setParent(currentNode);
                        if(!openSet.contains(neighbor)){openSet.add(neighbor);}
                        else{
                       	    openSet.remove(neighbor);
                            openSet.add(neighbor);
                        }
                    }                        
                }
            }
        }
        return null;
    }
    
    private ArrayList<int[]> retrace(Node4D startNode, Node4D endNode){
        ArrayList<int[]> path = new ArrayList<>();
        Node4D currentNode = endNode;
        while(currentNode != startNode){
        	int[] pos = currentNode.getPosition();
            path.add(new int[] {pos[0], pos[1], pos[2]});
            currentNode = currentNode.getParent();
        }
        Collections.reverse(path);
    	int[] pos = startNode.getPosition();
        path.add(0,new int[] {pos[0], pos[1], pos[2]});        
        return path;
    }    
}