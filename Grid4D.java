import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author Sara Azzouz <sazzouz@umich.edu>
 */

public class Grid4D {

	private int[] gridSize;
	private HashMap<String, Node4D> grid = new HashMap<>();

	public Grid4D(int width, int height, int nslice){
		gridSize = new int[] {width, height, nslice};
	}      

	public Node4D getNode(int xpos, int ypos, int zpos){
		String nodeKey = xpos+" "+ypos+" "+zpos;
		if(grid.containsKey(nodeKey))
			return grid.get(nodeKey);
		else if(isValid(xpos, ypos, zpos)) {
			grid.put(nodeKey, new Node4D(new int[]{xpos,ypos,zpos}));
			return grid.get(nodeKey);
		}
		return null;
	}

	public int getDistance(Node4D nodeA, Node4D nodeB){
		int[] posA = nodeA.getPosition();
		int[] posB = nodeB.getPosition();
		int dstX   = Math.abs(posA[0]-posB[0]);
		int dstY   = Math.abs(posA[1]-posB[1]);
		int dstZ   = Math.abs(posA[2]-posB[2]);
		int dstMin = Math.min(dstZ, Math.min(dstX, dstY));
		int dstMax = Math.max(dstZ, Math.max(dstX, dstY));
		int dstMid = dstX + dstY + dstZ - dstMin - dstMax;
		return (int)(10*((Math.sqrt(3)-Math.sqrt(2))*dstMin+(Math.sqrt(2)-1)*dstMid+dstMax));
	}

	public ArrayList<Node4D> getNeighbors(Node4D node){
		ArrayList<Node4D> neighbors = new ArrayList<>();
		for(int x = -1; x <= 1; ++x){
			for(int y = -1; y <= 1; ++y){
				for(int z = -1; z <= 1; ++z){
					if(x == 0 && y == 0 && z == 0) {continue;}
					int xpos = node.getPosition()[0]+x;
					int ypos = node.getPosition()[1]+y;
					int zpos = node.getPosition()[2]+z;
					if(isValid(xpos,ypos,zpos))
						neighbors.add(getNode(xpos,ypos,zpos));    			
				}
			}
		}
		return neighbors;
	}    

	public boolean isValid(int xpos, int ypos, int zpos){
		return xpos >= 0 && ypos >= 0 && zpos >= 0 && xpos < gridSize[0] && ypos < gridSize[1] && zpos < gridSize[2];
	}
} 