import java.awt.*;
import Jama.*; 
import java.util.ArrayList;
import java.util.HashMap;


public class MyPolygon extends Polygon {
	// For triangulation 
	public Polygon[] triangles;
	public int[] tempX = new int[3];
	public int[] tempY = new int[3];
	public int numTriangles; 
	public boolean triangulated = false; 

	// Mesh 
	public int[][] initialVertices; 
	public int[][] deformedVertices; 
	public int numVertices;
	public ArrayList<Integer> constrainedIdx; 
	public double[][][] triangleLocal; // Vertices in triangle local coordinates

	// public int[][] gComponentsSingle = new int[6][6];
	// public int[] constraintX = new int[3];
	// public int[] constraintY = new int[3]; 
	// public int numConstraints = 0; 
	// public int[][] gMatrix; 
	// public int[] v;
	// public Matrix gPrimeInvB;
	// public int[] newPos; 

	private void resizeTriangles(int num) {
		assert(tempX.length == tempY.length); 
		while(3*num > tempX.length) {
			// Double array capacity 
			int[] copyX = new int[tempX.length*2]; 
			int[] copyY = new int[tempX.length*2];
			// Copy data across 
			System.arraycopy(tempX, 0, copyX, 0, tempX.length);
			System.arraycopy(tempY, 0, copyY, 0, tempX.length);
			// Reset references 
			tempX = copyX; 
			tempY = copyY; 
		}
		assert(tempX.length == tempY.length);
	}

	/*
	private void resizeConstraints(int num) {
		assert(constraintX.length == constraintY.length);
		while(num > constraintX.length) {
			// Double array capacity
			int[] copyOldX = new int[constraintX.length*2]; 
			int[] copyOldY = new int[constraintX.length*2];
			// Copy data across 
			System.arraycopy(constraintX, 0, copyOldX, 0, tempX.length);
			System.arraycopy(constraintY, 0, copyOldY, 0, tempX.length);
			
			// Reset references 
			constraintX = copyOldX;
			constraintY = copyOldY;
		}
		assert(constraintX.length == constraintY.length);
	}
	*/

	public MyPolygon() {}
	public MyPolygon(int x[], int y[], int n) {}

	// Left 
	public boolean left(int ax, int ay, int bx, int by, int cx, int cy) {
		return (bx - ax)*(cy - ay) - (cx - ax)*(by - ay) < 0;
	}

	// Dot 
	public int DotProd(int[] x, int[] y) {
		return (x[0]*y[0]) + (x[1]*y[1]);
	}

	// Square length 
	public double squareLength(int[] vec) {
		return ((double)vec[0]*vec[0]) + ((double)vec[1]*vec[1]);
	}

	// O(n^4) implementation 
	public void triangulate() {
		numTriangles = 0; 

		// Compute z = x^2 + y^2
		int[] z = new int[npoints]; 
		for (int i=0; i<npoints; i++) {
			z[i] = (xpoints[i]*xpoints[i]) +
				(ypoints[i]*ypoints[i]);
		}
		// For each triple (i,j,k)
		int xn, yn, zn;
		double cx, cy;
		boolean flag; 
		for (int i=0; i<npoints-2; i++) {
			for (int j=i+1; j<npoints; j++) {
				for (int k=i+1; k<npoints; k++) {
					// Compute normal 
					xn = (ypoints[j]-ypoints[i])*(z[k]-z[i]) -
						(ypoints[k]-ypoints[i])*(z[j]-z[i]);
					yn = (xpoints[k]-xpoints[i])*(z[j]-z[i]) -
						(xpoints[j]-xpoints[i])*(z[k]-z[i]);
					zn = (xpoints[j]-xpoints[i])*(ypoints[k]-ypoints[i]) -
						(xpoints[k]-xpoints[i])*(ypoints[j]-ypoints[i]);
					// Examine faces on the bottom of the paraboloid
					if(zn < 0) {
						flag = true; 
						for(int m=0; m<npoints; m++) {
							flag = flag &&
								((xpoints[m]-xpoints[i])*xn +
								(ypoints[m]-ypoints[i])*yn +
								(z[m]-z[i])*zn <= 0);
						}
						if (flag) {
							// Exclude triangles entirely outside 
							cx = (xpoints[i] + xpoints[j] + xpoints[k]) / 3.0;
							cy = (ypoints[i] + ypoints[j] + ypoints[k]) / 3.0;
							if(contains(cx,cy)) {
								numTriangles++; 
								resizeTriangles(numTriangles);
								tempX[3*(numTriangles-1)] = xpoints[i];
								tempX[3*(numTriangles-1) + 1] = xpoints[j];
								tempX[3*(numTriangles-1) + 2] = xpoints[k];
								tempY[3*(numTriangles-1)] = ypoints[i];
								tempY[3*(numTriangles-1) + 1] = ypoints[j];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
							}
						}
					}
				}
			}
		}
		triangles = new Polygon[numTriangles];
		int[] newPolyX = new int[3]; 
		int[] newPolyY = new int[3];
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++) {
				newPolyX[j] = tempX[3*i + j];
				newPolyY[j] = tempY[3*i + j]; 
			}
			triangles[i] = new Polygon(newPolyX, newPolyY, 3);
		}
		triangulated = true; 
	}

	// Initialize from Mesh
	public void initializeMesh() {
		constrainedIdx = new ArrayList<Integer>();
		constrainedIdx.clear();

		// Copy vertices 
		ArrayList<int[]> verts = new ArrayList<int[]>(); 
		boolean inVerts; 
		numVertices = 0; 
		for(int i=0; i<numTriangles; i++) { 
			Polygon triangle = triangles[i]; 
			for(int j=0; j<3; j++) {
				inVerts = false;
				for(int k=0; k<numVertices; k++) {
					if(triangle.xpoints[j] == verts.get(k)[0] && 
						triangle.ypoints[j] == verts.get(k)[1]) {
						inVerts = true; 
						break;
					}
				}
				if(!inVerts) {
					verts.add(new int[]{triangle.xpoints[j], triangle.ypoints[j]});
					numVertices++;
				}
			}
		}
		initialVertices = new int[numVertices][2];
		deformedVertices = new int[numVertices][2];
		for(int i=0; i<numVertices; i++) { 
			int[] aVertex = verts.get(i);
			initialVertices[i][0] = aVertex[0];
			initialVertices[i][1] = aVertex[1];
			deformedVertices[i][0] = aVertex[0];
			deformedVertices[i][1] = aVertex[1];
		}

		// Triangle in local coodinates 
		triangleLocal = new double[numTriangles][3][2]; 
		int n0, n1, n2; 
		int[] v0 = new int[2], v1 = new int[2], v2 = new int[2];
		int[] v01 = new int[2], v01Rot90 = new int[2];
		int[] v02 = new int[2];
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++) {
				n0 = j; 
				n1 = (j+1)%3;
				n2 = (j+2)%3; 

				v0[0] = triangles[i].xpoints[0];
				v0[1] = triangles[i].ypoints[0];
				v1[0] = triangles[i].xpoints[1];
				v1[1] = triangles[i].ypoints[1];
				v2[0] = triangles[i].xpoints[2];
				v2[1] = triangles[i].ypoints[2];

				v01[0] = v1[0] - v0[0];
				v01[1] = v1[1] - v0[0];
				v01Rot90[0] = v01[1];
				v01Rot90[1] = -v01[0];

				v02[0] = v2[0] - v0[0];
				v02[1] = v2[1] - v0[1];

				double fx = (double)DotProd(v02, v01) / squareLength(v01);
				double fy = (double)DotProd(v02, v01Rot90) / squareLength(v01Rot90);
				triangleLocal[i][j][0] = fx; 
				triangleLocal[i][j][1] = fy; 

				// Sanity check 
				System.out.println(v0[0] + fx*v01[0] + fy*v01Rot90[0] - v2[0]);
				System.out.println(v0[1] + fx*v01[1] + fy*v01Rot90[1] - v2[1]);
			}
		}
	}


	public int[] addConstraint(int x, int y) {
		int distToPoint;
		int[] aConstraint; 
		for(int i=0; i<constrainedIdx.size(); i++) {
			aConstraint = deformedVertices[constrainedIdx.get(i)]; 
			distToPoint = (aConstraint[0]-x)*(aConstraint[0]-x) +
				(aConstraint[1]-y)*(aConstraint[1]-y);
			if(distToPoint < 25) {
				return new int[]{0, 0, 0}; 
			}
		}
		for(int i=0; i<numVertices; i++) { 
			distToPoint = (deformedVertices[i][0]-x)*(deformedVertices[i][0]-x) +
				(deformedVertices[i][1]-y)*(deformedVertices[i][1]-y);
			if(distToPoint < 25) {
				constrainedIdx.add(i);
				return new int[]{1, deformedVertices[i][0], deformedVertices[i][1]}; 
			}
		}

		return new int[]{0, 0, 0};
	}

	public void calcStepOneMatrix() { 
		double[][] gMatrix = new double[2*numVertices][2*numVertices];
		// Zero 
		for(int i=0; i<2*numVertices; i++) {
			for(int j=0; j<2*numVertices; j++) {
				gMatrix[i][j] = 0.0;
			}
		}

		// Number of constraints 
		int numConstraints = constrainedIdx.size(); 
		int numFreeVars = numVertices - numConstraints;

		// Figure out vertex ordering 
		HashMap<XYKey, Integer> vertMap = new HashMap<XYKey, Integer>();

		boolean isConstraint; 
		int nRow = 0; 
		for(int i=0; i<numVertices; i++) {
			isConstraint = false; 
			for(int j=0; j<numConstraints; j++) {
				if(constrainedIdx.get(j) == i) {
					isConstraint = true; 
					break;
				}
			}
			if(!isConstraint) {
				//System.out.println(i);
				//System.out.println(initialVertices[i][0] + ", " + initialVertices[i][1]);
				vertMap.put(new XYKey(initialVertices[i][0],initialVertices[i][1]), nRow);
				nRow++; 
			} 
		}
		for(int i=0; i<numConstraints; i++) {
			//System.out.println(constrainedIdx.get(i));
			//System.out.println(initialVertices[constrainedIdx.get(i)][0] +
			//	", " + initialVertices[constrainedIdx.get(i)][1]);
			vertMap.put(new XYKey(initialVertices[constrainedIdx.get(i)][0], 
				initialVertices[constrainedIdx.get(i)][1]), nRow);
			nRow++; 
		}
		//System.out.println(vertMap.size());
		//System.out.println(initialVertices.length);
				

		// Build vector for computing the G matrix 
		int[] vVector = new int[2*numVertices];
		for(int i=0; i<numVertices; i++) {
			isConstraint = false; 
			for(int j=0; j<constrainedIdx.size(); j++) {
				if(i == constrainedIdx.get(j)) {
					isConstraint = true; 
					break; 
				}
			} 
			if(!isConstraint) {
				int rowInVec = vertMap.get(new XYKey(initialVertices[i][0], initialVertices[i][1]));
				vVector[2*rowInVec] = initialVertices[i][0];
				vVector[2*rowInVec + 1] = initialVertices[i][1]; 
			}
		}
		for(int i=0; i<numConstraints; i++) {
			int rowInVec = vertMap.get(
				new XYKey(initialVertices[constrainedIdx.get(i)][0], initialVertices[constrainedIdx.get(i)][1]));
			vVector[2*rowInVec] = initialVertices[constrainedIdx.get(i)][0];
			vVector[2*rowInVec + 1] = initialVertices[constrainedIdx.get(i)][1]; 
		}

		// Build the G Matrix 
		int n0x, n0y, n1x, n1y, n2x, n2y; 
		double x, y; 
		for(int i=0; i<numTriangles; i++) {
			Polygon triangle = triangles[i]; 
			for(int j=0; j<3; j++) {
				n0x = 2 * vertMap.get(new XYKey(triangle.xpoints[j], triangle.ypoints[j]));
				n0y = n0x + 1; 
				n1x = 2 * vertMap.get(new XYKey(triangle.xpoints[(j+1)%3], triangle.ypoints[(j+1)%3]));
				n1y = n1x + 1; 
				n2x = 2 * vertMap.get(new XYKey(triangle.xpoints[(j+2)%3], triangle.ypoints[(j+2)%3]));
				n2y = n2x + 1; 
				x = triangleLocal[i][j][0];
				y = triangleLocal[i][j][1];
				System.out.println(n0x);
				System.out.println(n1x);
				System.out.println(n2x);

				// Sanity check 
				int[] v0 = new int[]{vVector[n0x], vVector[n0y]};
				int[] v1 = new int[]{vVector[n1x], vVector[n1y]};
				int[] v2 = new int[]{vVector[n2x], vVector[n2y]};
				int[] v01 = new int[]{v1[0] - v0[0], v1[1] - v0[1]};
				int[] v01Perp = new int[]{v01[1], -v01[0]};
				
				
			}
		}

	}
	

	/*
	public int[] setConstraintActive(int x, int y) {
		int distToPoint; 
		for(int i=0; i<numConstraints; i++) {
			distToPoint = (constraintX[i] - x)*(constraintX[i] - x) + (constraintY[i] - y)*(constraintY[i] - y); 
			if(distToPoint < 25) { 
				return new int[]{1, i};
			} 
		}
		return new int[]{0, 0};
	}
	

	// Precompute components for scale-free manipulation
	// Step One in Igarashi et. al's paper
	// Helper method for calculating 
	// matrix components for a single triangle
	// @params: v0x, v0y, v1x, v1y, v2x, v2y 
	// Output: integer array of matrix components 
	public void calcComponentsSingle(int v0x, int v0y, int v1x, int v1y, int v2x, int v2y) {
		// find x01, y01
		int x01 = (v2x-v0x)*(v1x-v0x) + (v2y-v0y)*(v1y-v0y);
		int y01 = (v2x-v0x)*(v1y-v0y) + (v2y-v0y)*(v0x-v1x);

		// Populate gComponentsSingle[] array 
		for(int i=0; i<gComponentsSingle.length; i++){
			for(int j=0; j<gComponentsSingle[0].length; j++) {
				gComponentsSingle[i][j] = 0; 
			}
		}
		gComponentsSingle[0][0] = x01*x01 - 2*x01 + y01*y01 + 1; 
		gComponentsSingle[2][0] = -2*x01*x01 + 2*x01 - 2*y01*y01; 
		gComponentsSingle[3][0] = 2*y01; 
		gComponentsSingle[4][0] = 2*x01 - 2; 
		gComponentsSingle[5][0] = -2*y01; 
		gComponentsSingle[1][1] = x01*x01 - 2*x01 + y01*y01 + 1; 
		gComponentsSingle[1][2] = -2*y01; 
		gComponentsSingle[1][3] = -2*x01*x01 + 2*x01 - 2*y01*y01; 
		gComponentsSingle[1][4] = 2*y01; 
		gComponentsSingle[1][5] = 2*x01 - 2; 
		gComponentsSingle[2][2] = x01*x01 + y01*y01; 
		gComponentsSingle[2][4] = -2*x01;
		gComponentsSingle[2][5] = 2*y01;
		gComponentsSingle[3][3] = x01*x01 + y01*y01; 
		gComponentsSingle[3][4] = -2*y01; 
		gComponentsSingle[3][5] = -2*x01; 
		gComponentsSingle[4][4] = 1; 
		gComponentsSingle[5][5] = 1;
	}
	// Helper method for populating the GMatrix 
	public void buildG(int idx0, int idx1, int idx2) {
		int[] actual = new int[]{idx0, idx0+1, idx1, idx1+1, idx2, idx2+1}; 
		for(int i=0; i < gComponentsSingle.length; i++) {
			for(int j=0; j <gComponentsSingle[0].length; j++) {
				gMatrix[actual[i]][actual[j]] += gComponentsSingle[i][j];
			}
		}
	}

	public void calcGMatrix() { 
		// Iterate through triangles 
		// Find unique verticies 
		ArrayList<Integer> vListX = new ArrayList<Integer>(); 
		ArrayList<Integer> vListY = new ArrayList<Integer>();
		int x, y; 
		boolean inV;
		for(int i=0; i<numTriangles; i++) { 
			for(int j=0; j<triangles[i].npoints; j++) {
				x = triangles[i].xpoints[j];
				y = triangles[i].ypoints[j];
				inV = false; 
				for(int k=0; k<vListX.size(); k++) {
					if(x == vListX.get(k) && y == vListY.get(k)){
						inV = true; 
						break;
					}
				}
				if(!inV) {
					vListX.add(x);
					vListY.add(y);
				}
			}
		}
		v = new int[vListX.size() * 2]; 
		for(int i=0; i<vListX.size(); i++){
			v[2*i] = vListX.get(i); 
			v[(2*i) + 1] = vListY.get(i);
		}

		gMatrix = new int[v.length][v.length];
		for(int i=0; i<v.length; i++) {
			for(int j=0; j<v.length; j++) {
				gMatrix[i][j] = 0;
			}
		}
		int[] vIdx = new int[3];
		for(int i=0; i<numTriangles; i++) { 
			// Get triangle vertex indices in v 
			for(int j=0; j<triangles[i].npoints; j++) {
				for(int k=0; k<v.length; k+=2) {
					if(triangles[i].xpoints[j] == v[k] && triangles[i].ypoints[j] == v[k+1]) {
						vIdx[j] = k;
						break;
					}
				}
			}
			// Calculate associated G components 
			calcComponentsSingle(triangles[i].xpoints[0], triangles[i].ypoints[0], 
				triangles[i].xpoints[1], triangles[i].ypoints[1], 
				triangles[i].xpoints[2], triangles[i].ypoints[2]);
			buildG(vIdx[0], vIdx[1], vIdx[2]);

			calcComponentsSingle(triangles[i].xpoints[1], triangles[i].ypoints[1], 
				triangles[i].xpoints[2], triangles[i].ypoints[2], 
				triangles[i].xpoints[0], triangles[i].ypoints[0]);
			buildG(vIdx[1], vIdx[2], vIdx[0]);

			calcComponentsSingle(triangles[i].xpoints[2], triangles[i].ypoints[2], 
				triangles[i].xpoints[0], triangles[i].ypoints[0], 
				triangles[i].xpoints[1], triangles[i].ypoints[1]);
			buildG(vIdx[2], vIdx[0], vIdx[1]);
		}
	}

	// Function for obtaining G' and B in scale-free construction 
	public void getGPrimeInvB() {
		ArrayList<Integer> constraintIdx = new ArrayList<Integer>(); 
		ArrayList<Integer> notConstraintIdx = new ArrayList<Integer>(); 
		boolean isConstraint; 
		for(int i=0; i<v.length; i+=2) { 
			isConstraint = false; 
			for(int j=0; j<numConstraints; j++) { 
				if(v[i] == constraintX[j] && v[i+1] == constraintY[j]) {
					isConstraint = true; 
					break;
				} 
			}
			if(isConstraint) {
				constraintIdx.add(i);
			} else {
				notConstraintIdx.add(i); 
			}
		}

		int[][] gReorderRow = new int[gMatrix.length][gMatrix.length]; 
		for(int i=0; i<numConstraints; i++) {
			for(int j=0; j<gMatrix.length; j++) {
				gReorderRow[gMatrix.length - 2*numConstraints + i][j] = gMatrix[constraintIdx.get(i)][j];
				gReorderRow[gMatrix.length - 2*numConstraints + i + 1][j] = gMatrix[constraintIdx.get(i) + 1][j];
			}
		}
		int numFreeVariables = (gMatrix.length / 2) - numConstraints;
		for(int i=0; i<numFreeVariables; i++) {
			for(int j=0; j<gMatrix.length; j++) {
				gReorderRow[i][j] = gMatrix[notConstraintIdx.get(i)][j];
				gReorderRow[i + 1][j] = gMatrix[notConstraintIdx.get(i) + 1][j];
			}
		}
		
		int[][] gReorder = new int[gMatrix.length][gMatrix.length];
		for(int i=0; i<gMatrix.length; i++) {
			for(int j=0; j<numConstraints; j++) {
				gReorder[i][gMatrix.length - 2*numConstraints + j] = gReorderRow[i][constraintIdx.get(j)];
				gReorder[i][gMatrix.length - 2*numConstraints + j + 1] = gReorderRow[i][constraintIdx.get(j) + 1];
			}
		}
		for(int i=0; i<gMatrix.length; i++) {
			for(int j=0; j<numFreeVariables; j++) {
				gReorder[i][j] = gReorderRow[i][notConstraintIdx.get(j)];
				gReorder[i][j + 1] = gReorderRow[i][notConstraintIdx.get(j) + 1];
			}
		}

		double[][] g00 = new double[gMatrix.length - 2*numConstraints][gMatrix.length - 2*numConstraints];
		double[][] g10 = new double[2*numConstraints][gMatrix.length - 2*numConstraints];
		double[][] g01 = new double[gMatrix.length - 2*numConstraints][2*numConstraints];

		for(int i=0; i<gMatrix.length - 2*numConstraints; i++) {
			for(int j=0; j<gMatrix.length - 2*numConstraints; j++) {
				g00[i][j] = (double)gMatrix[i][j];
			}
		}
		for(int i=0; i<2*numConstraints; i++) {
			for(int j=0; j<gMatrix.length - 2*numConstraints; j++) {
				g10[i][j] = (double)gMatrix[gMatrix.length - 2*numConstraints + i][j];
			}
		}
		for(int i=0; i<gMatrix.length - 2*numConstraints; i++) {
			for(int j=0; j<2*numConstraints; j++) {
				g01[i][j] = (double)gMatrix[i][gMatrix.length - 2*numConstraints + j];
			}
		}

		Matrix g00Matrix = new Matrix(g00);
		Matrix g10Matrix = new Matrix(g10);
		Matrix g01Matrix = new Matrix(g01);
		
		Matrix gPrime = g00Matrix.plus(g00Matrix.transpose());
		Matrix b = g01Matrix.plus(g10Matrix.transpose());
		b = b.times(-1);
		Matrix gPrimeInv = gPrime.inverse(); 
		gPrimeInvB = gPrimeInv.times(b);
	}
	
	// Calculating new coordinates 
	public void shapeManipulate() {
		//System.out.println(gPrimeInvB.getRowDimension() + ", " + gPrimeInvB.getColumnDimension());
		double[][] constraints = new double[2*numConstraints][1];
		for(int i=0; i<numConstraints; i++) { 
			constraints[2*i][0] = constraintX[i];
			constraints[2*i + 1][0] = constraintY[i];
		}
		Matrix qMatrix = new Matrix(constraints); 
		//System.out.println(qMatrix.getRowDimension() + ", " + qMatrix.getColumnDimension());
		Matrix manipulateNewPos = gPrimeInvB.times(qMatrix);
		double[][] newPosDouble = manipulateNewPos.getArrayCopy();
		//System.out.println(newPosDouble[0].length);
		newPos = new int[newPosDouble.length]; 
		System.out.println("v length: " + v.length);
		System.out.println("New Length: " + newPosDouble.length);

		for(int i=0; i<newPosDouble.length; i+=2) {
			System.out.println("Old: " + v[i] + ", " + v[i+1]);
			System.out.println("New: " + (int)newPosDouble[i][0] + ", " + (int)newPosDouble[i+1][0]);
		}
	}
	*/
}