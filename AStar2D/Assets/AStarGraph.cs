using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;

/*
SOURCES:
https://docs.unity3d.com/ScriptReference/NavMesh.html
https://msdn.microsoft.com/en-us/library/aa289152(v=vs.71).aspx
http://blog.two-cats.com/2014/06/a-star-example/
https://en.wikipedia.org/wiki/A*_search_algorithm
https://www.youtube.com/watch?v=mZfyt03LDH4
https://www.youtube.com/watch?v=3Dw5d7PlcTM
*/

public class Heap<T> where T : IHeapItem<T>
{
    T[] items;
    int heapCount;

    public Heap(int maxSize)
    {
        items = new T[maxSize];
    }

    public void Add(T item)
    {
        item.HeapIndex = heapCount;
        items[heapCount] = item;
        bubbleUp(item);
        heapCount++;
    }

    public bool Contains(T item)
    {
        return Equals(items[item.HeapIndex], item);
    }

    public T removeFirst()
    {
        T firstItem = items[0];
        heapCount--;
        items[0] = items[heapCount];
        items[0].HeapIndex = 0;
        bubbleDown(items[0]);
        return firstItem;
    }

    void bubbleDown(T item)
    {
        while (true)
        {
            int childLeft = item.HeapIndex * 2 + 1;
            int childRight = item.HeapIndex * 2 + 2;
            int swapIndex = 0;

            if (childLeft < heapCount)
            {
                swapIndex = childLeft;

                if (childRight < heapCount)
                {
                    if (items[childLeft].CompareTo(items[childRight]) < 0)
                    {
                        swapIndex = childRight;
                    }
                }

                if (item.CompareTo(items[swapIndex]) < 0)
                    Swap(item, items[swapIndex]);
                else
                    return;
            }
            else
                return;
        }
    }

    void bubbleUp(T item)
    {
        int pidx = (item.HeapIndex - 1) / 2;

        while (true)
        {
            T parent = items[pidx];
            if (item.CompareTo(parent) > 0)
                Swap(item, parent);
            else
                break;
            pidx = (item.HeapIndex - 1) / 2;
        }
    }

    public int Count
    {
        get
        {
            return heapCount;
        }
    }

    void Swap(T itemOne, T itemTwo)
    {
        items[itemOne.HeapIndex] = itemTwo;
        items[itemTwo.HeapIndex] = itemOne;
        int itemOneIndex = itemOne.HeapIndex;
        itemOne.HeapIndex = itemTwo.HeapIndex;
        itemTwo.HeapIndex = itemOne.HeapIndex;
    }

}

public interface IHeapItem<T> : IComparable<T>
{
    int HeapIndex
    {
        get;
        set;
    }
}


public class Node : IHeapItem<Node>
{
    object data;
    public string Key;
    //public AdjacencyList neighbors;
    public List<Edge> neighbors;
    public Vector3 worldPos;

    public Node parent;
    public float GCost;
    public float Heuristic;
    public float FCost;
    int heapIndex;

    public Node(string key, object data, Vector3 worldPos)
    {
        this.data = data;
        this.Key = key;

        //this.neighbors = new AdjacencyList();
        this.neighbors = new List<Edge>();

        //this.parent = null;
        this.FCost = 0.0f;
        this.Heuristic = 0.0f;
        this.GCost = 0.0f;
        this.worldPos = worldPos;
    }

    protected internal virtual void AddDirected(Node n) { AddDirected(new Edge(n, 0)); }
    protected internal virtual void AddDirected(Node n, int cost) { AddDirected(new Edge(n, cost)); }
    protected internal virtual void AddDirected(Edge e)
    {
        neighbors.Add(e);
        //Debug.Log("Added "+e.toString()+" to "+this.Key);
    }

    public int HeapIndex
    {
        get { return heapIndex; }
        set { heapIndex = value; }
    }

    public int CompareTo(Node ntoComp)
    {
        int comp = FCost.CompareTo(ntoComp.FCost);
        if (comp == 0)
        {
            comp = Heuristic.CompareTo(ntoComp.Heuristic);
        }
        return -comp;
    }

    public string toString()
    {
        return "Node" + this.Key + ":" + data.ToString();
    }
}


public class SearchNode
{
    private Node reference;
    private double H;
    private double G;

    public SearchNode()
    {

    }
}


public class SearchTree
{
    private Dictionary<SearchNode,SearchNode> searchMap;

    public SearchTree()
    {
        searchMap = new Dictionary<SearchNode, SearchNode>();
    }
}


public class Edge
{
    protected int cost;
    protected Node neighbor;

    public Edge(Node neighbor) : this(neighbor, 0) { }

    public Edge(Node neigh, int cos)
    {
        cost = cos;
        neighbor = neigh;
    }

    public Node getNeighbor() { return this.neighbor; }
    public int getCost() { return this.cost; }

    public string toString()
    {
        string costString = cost.ToString();
        string neigh = neighbor.toString();
        return neigh+" : "+costString;
    }
}


public class AdjacencyList : CollectionBase
{
    protected internal virtual void Add(Edge e)
    {
        base.InnerList.Add(e);
    }

    public virtual Edge this[int index]
    {
        get { return (Edge)base.InnerList[index]; }
        set { base.InnerList[index] = value; }
    }
}


public class NodeList// : IEnumerable
{
    public Hashtable data;
    private int size;

    public NodeList()
    {
        data = new Hashtable();
        size = 0;
    }

    public virtual void Add(Node n)
    {
        data.Add(n.Key, n);
        this.size++;
    }

    public virtual void Remove(Node n)
    {
        data.Remove(n.Key);
        this.size--;
    }

    public int getSize()
    { return this.size; }

    public virtual bool ContainsKey(string key)
    {
        return data.ContainsKey(key);
    }

    public virtual void Clear()
    {
        data.Clear();
        this.size = 0;
    }

    // Properties...
    public virtual Node this[string key]
    {
        get
        {
            return (Node)data[key];
        }
    }
}


public class Graph
{
    private NodeList nodeList;

    public Graph()
    {
        this.nodeList = new NodeList();
    }

    public virtual Node AddNode(string s, object o, Vector3 worldPos)
    {
        if (!nodeList.ContainsKey(s))
        {
            Node n = new Node(s, o, worldPos);
            nodeList.Add(n);
            //Debug.Log("IT'S POS " + n.worldPos);
            return n;
        }
        else
        {
            throw new System.Exception("There is already a node of that type: "+s);
        }
    }

    public void AddNode(Node n)
    {
        if (!nodeList.ContainsKey(n.Key))
        {
            nodeList.Add(n);
        }
        else
        {
            throw new System.Exception("There is already a node of that type: "+n.Key);
        }
    }

    public virtual void AddDirectedEdge(string uKey, string vKey)
    {
        AddDirectedEdge(uKey, vKey, 0);
    }

    public virtual void AddDirectedEdge(string uKey, string vKey, int cost)
    {
        // get references to uKey and vKey
        if (nodeList.ContainsKey(uKey) && nodeList.ContainsKey(vKey))
            AddDirectedEdge(nodeList[uKey], nodeList[vKey], cost);
        else
            throw new System.Exception("Error adding directed edge "+uKey+" : "+vKey);
    }

    public virtual void AddDirectedEdge(Node u, Node v)
    {
        AddDirectedEdge(u, v, 0);
    }

    public virtual void AddDirectedEdge(Node u, Node v, int cost)
    {
        // Make sure u and v are Nodes in this graph
        if (nodeList.ContainsKey(u.Key) && nodeList.ContainsKey(v.Key))
            // add an edge from u -> v
            u.AddDirected(v, cost);
        else
            // one or both of the nodes were not found in the graph
            throw new System.Exception("Error adding directed edge");
    }


    public virtual void AddUndirectedEdge(string uKey, string vKey)
    {
        AddUndirectedEdge(uKey, vKey, 0);
    }

    public virtual void AddUndirectedEdge(string uKey, string vKey, int cost)
    {
        // get references to uKey and vKey
        if (nodeList.ContainsKey(uKey) && nodeList.ContainsKey(vKey))
            AddUndirectedEdge(nodeList[uKey], nodeList[vKey], cost);
        else
            throw new System.Exception("Error adding undirected edge");
    }

    public virtual void AddUndirectedEdge(Node u, Node v)
    {
        AddUndirectedEdge(u, v, 0);
    }

    public virtual void AddUndirectedEdge(Node u, Node v, int cost)
    {
        // Make sure u and v are Nodes in this graph
        if (nodeList.ContainsKey(u.Key) && nodeList.ContainsKey(v.Key))
        {
            // Add an edge from u -> v and from v -> u
            u.AddDirected(v, cost);
            v.AddDirected(u, cost);
        }
        else
            // one or both of the nodes were not found in the graph
            throw new System.Exception("Error adding undirected edge");
    }


    public virtual bool Contains(Node n)
    {
        return Contains(n.Key);
    }

    public virtual bool Contains(string key)
    {
        return nodeList.ContainsKey(key);
    }

    public NodeList returnNodes() { return this.nodeList; }
}


public class AStar
{
    protected double Heuristic;
    protected double G;
    protected double F;

    protected Graph graph;

    public List<Node> path;
    private bool useHeap;

    public AStar(Graph graph, bool useHeap)
    {
        this.graph = graph;
        this.path = new List<Node>();
        this.useHeap = useHeap;
    }

    protected void calculateF()
    { this.F = this.G + this.Heuristic; }

    public void findPath(Node startNode, Node endNode)
    {
        Debug.Log("Finding Path.");

        List<Node> openSet = new List<Node>();

        Heap<Node> openSetHeap = new Heap<Node>(2000);

        HashSet<Node> closedSet = new HashSet<Node>();
        openSet.Add(startNode);
        openSetHeap.Add(startNode);

        if (!useHeap)
        {
            while (openSet.Count > 0)
            {
                Node currentNode = openSet[0];

                for (int i = 0; i < openSet.Count; ++i)
                {
                    if (openSet[i].FCost < currentNode.FCost || openSet[i].FCost == currentNode.FCost)
                    {
                        if (openSet[i].Heuristic < currentNode.Heuristic)
                        {
                            currentNode = openSet[i];
                            Debug.Log("Current Node = " + openSet[i].toString());
                        }
                    }
                }

                openSet.Remove(currentNode);
                closedSet.Add(currentNode);

                if (currentNode == endNode)// End Case - We're there
                {
                    RetracePath(startNode, endNode);
                    Debug.Log("Finished.");
                    return;
                }

                foreach (Edge e in currentNode.neighbors)
                {
                    Node neighbor = e.getNeighbor();

                    //Debug.Log("Processing neighbor "+ neighbor.toString());
                    if (closedSet.Contains(neighbor))
                    {
                        continue;
                    }

                    float MovementCost = currentNode.GCost + getDistance(currentNode, neighbor);
                    if (MovementCost < neighbor.GCost || !openSet.Contains(neighbor))
                    {
                        neighbor.GCost = MovementCost;
                        neighbor.Heuristic = getDistance(neighbor, endNode);
                        neighbor.parent = currentNode;

                        if (!openSet.Contains(neighbor))
                        {
                            openSet.Add(neighbor);
                            //Debug.Log("Adding Neighbor. " + neighbor.Key);
                        }
                    }
                }
            }
        }
        else
        {
            while (openSetHeap.Count > 0)
            {
                Node currentNode = openSetHeap.removeFirst();
                Debug.Log("Current Node: "+currentNode.Key);
                closedSet.Add(currentNode);

                if (currentNode == endNode)// End Case - We're there
                {
                    RetracePath(startNode, endNode);
                    Debug.Log("Finished.");
                    return;
                }

                foreach (Edge e in currentNode.neighbors)
                {
                    Node neighbor = e.getNeighbor();

                    if (closedSet.Contains(neighbor))
                    {
                        continue;
                    }

                    float MovementCost = currentNode.GCost + getDistance(currentNode, neighbor);
                    if (MovementCost < neighbor.GCost || !openSet.Contains(neighbor))
                    {
                        neighbor.GCost = MovementCost;
                        neighbor.Heuristic = getDistance(neighbor, endNode);
                        neighbor.parent = currentNode;

                        if (!openSetHeap.Contains(neighbor))
                        {
                            openSetHeap.Add(neighbor);
                            //Debug.Log("Adding Neighbor. " + neighbor.Key);
                        }
                    }
                }
            }
        }
    }

    private void RetracePath(Node startNode, Node endNode)
    {
        List<Node> pth = new List<Node>();
        Node currentNode = endNode;

        while (currentNode != startNode)
        {
            pth.Add(currentNode);
            currentNode = currentNode.parent;
        }
        pth.Reverse();
        this.path = pth;
        Debug.Log("Retraced Path");
    }

    private float getDistance(Node one, Node two)
    {
        float fone = Mathf.Abs(one.worldPos.x - two.worldPos.x); 
        float ftwo = Mathf.Abs(one.worldPos.z - two.worldPos.z);

        if (fone > ftwo)
            return 14* ftwo + 10*(fone - ftwo);
        return 14 * fone + 10 * (ftwo - fone);
    }
}


public class AStarGraph : MonoBehaviour
{
    private Vector3[] verts;
    private Graph g;

    private List<Vector3> nodePositions;
    private List<Node> nodeContainer;

    public bool showPaths;
    public bool useHeap;
    public Transform cylinder;
    public Transform endGuy;

    private int startNodeIndex = 30;
    private int endNodeIndex = 55;

    private AStar astar;
    private NavMeshTriangulation nmtri;

    protected Vector3 getCentroid(Vector3 vertex1, Vector3 vertex2, Vector3 vertex3)
    {
        float sumOne = vertex1.x + vertex2.x + vertex3.x;
        float sumTwo = vertex1.y + vertex2.y + vertex3.y;
        float sumThree = vertex1.z + vertex2.z + vertex3.z;

        return new Vector3(sumOne/3, sumTwo/3, sumThree/3);
    }

    // Use this for initialization
    void Start()
    {
        nodePositions = new List<Vector3>();
        nodeContainer = new List<Node>();

        verts = NavMesh.CalculateTriangulation().vertices;
        int[] indicies = NavMesh.CalculateTriangulation().indices;
        nmtri = NavMesh.CalculateTriangulation();

        g = new Graph();

        //Add Nodes
        Vector3 ctoid = new Vector3();
        for (int k = 0; k < nmtri.vertices.Length - 2; ++k)
        {
            ctoid = getCentroid(nmtri.vertices[k], nmtri.vertices[k + 1], nmtri.vertices[k + 2]);
            nodePositions.Add(ctoid);
            Node reff = g.AddNode("Node: " + k.ToString(), k, ctoid);
            nodeContainer.Add(reff);
        }

        Vector3 toid = new Vector3();
        int jj = verts.Length;
        toid = getCentroid(nmtri.vertices[jj - 3], nmtri.vertices[jj - 2], nmtri.vertices[jj - 1]);
        nodePositions.Add(toid);
        Node refff = g.AddNode("Node: " + (jj - 2).ToString(), jj-2, toid);
        nodeContainer.Add(refff);

        Debug.Log("Added Nodes to Graph.\n");
        Debug.Log(verts.Length);
        //Add edges
        
        /*
        for (int t = 0; t < nodeContainer.Count - 15; ++t)
        {
            for (int tt = 0; tt < nodeContainer.Count - (10); ++tt)
            {
                Node n = nodeContainer[t];
                Node no = nodeContainer[tt];
                if (n != no)
                {
                    double dist = Mathf.Sqrt(Mathf.Pow((n.worldPos.x - no.worldPos.x), 2) + Mathf.Pow((n.worldPos.z - no.worldPos.z), 2));
                    Debug.Log("DISTANCE: " + dist.ToString());
                    if (dist <= 45.0)
                    {
                        g.AddDirectedEdge("Node: "+t.ToString(), "Node: "+tt.ToString());
                        //Debug.Log("Added Edge "+ "Node: " + t.ToString()+" "+ "Node: " + tt.ToString());
                    }
                }
            }
        }
        */

        
        /**/
        for (int j = 0; j < nodeContainer.Count - 2; ++j)
        {
            g.AddDirectedEdge("Node: " + (j).ToString(), "Node: " + (j + 1).ToString());
        }
        g.AddDirectedEdge("Node: 47", "Node: 0");
        g.AddDirectedEdge("Node: 43", "Node: 0");
        /**/
        

        Debug.Log("Added Edges to Graph.\n");

        // Begin Astar Pathfinding
        astar = new AStar(g, useHeap);

        if (this.nodeContainer.Count == 0)
            throw new System.Exception("Problem getting node container info.?");

        cylinder.position = this.nodeContainer[startNodeIndex].worldPos;
        endGuy.position = this.nodeContainer[endNodeIndex].worldPos;

        astar.findPath(this.nodeContainer[startNodeIndex], this.nodeContainer[endNodeIndex]);
        Debug.Log("endNodeIndex" + endNodeIndex.ToString());
    }

    void OnDrawGizmos()
    {
        if (!Application.isPlaying)
            return;   
        Gizmos.color = Color.red;
        if (showPaths)
        {
            for (int i = 0; i < nodeContainer.Count; ++i)
            {
                Gizmos.DrawWireCube(nodeContainer[i].worldPos, new Vector3(0.5f, 0.5f, 0.5f));
            }
        }
        
        if (showPaths && astar.path.Count > 0)
        {
            Gizmos.color = Color.black;
            foreach (Node n in astar.path)
            {
                Gizmos.DrawCube(n.worldPos, new Vector3(1.0f, 1.0f, 1.0f));
            }
            //Gizmos.DrawCube(endGuy.position, new Vector3(0.25f, 0.25f, 0.25f));
        }
        
    }

    void redoAstar(bool rightArrow)
    {
        if (astar.path.Count > 0)
            astar.path.Clear();

        if (rightArrow)
        {
            if (endNodeIndex < this.nodeContainer.Count)
                endNodeIndex++;
        }
        else
        {
            if (endNodeIndex > 0)
                endNodeIndex--;
        }

        if (this.nodeContainer.Count == 0)
            throw new System.Exception("Problem getting node container info.?");

        cylinder.position = this.nodeContainer[startNodeIndex].worldPos;
        endGuy.position = this.nodeContainer[endNodeIndex].worldPos;

        astar.findPath(this.nodeContainer[startNodeIndex], this.nodeContainer[endNodeIndex]);
        Debug.Log("endNodeIndex" + endNodeIndex.ToString());
    }

    // Update is called once per frame
    void Update ()
    {
        /**/
        for (int i = 0; i < nodeContainer.Count - 1; ++i)
        {
            Debug.DrawLine(nodeContainer[i].worldPos, nodeContainer[i + 1].worldPos, Color.cyan);
        }
        /**/

        for (int i = 0; i < astar.path.Count - 1; ++i)
        {
            Debug.DrawLine(astar.path[i].worldPos, astar.path[i + 1].worldPos, Color.magenta);
        }
        /**/

        /*for (int i = 0; i < nmtri.indices.Length-1; ++i)
        {
            Debug.DrawLine(nmtri.vertices[nmtri.indices[i]], nmtri.vertices[nmtri.indices[i + 1]], Color.magenta);
        }*/

        //Draw Neighbors
        /*
        foreach (Edge n in nodeContainer[0].neighbors)
        {
            Node neigh = n.getNeighbor();
            Debug.DrawLine(nodeContainer[0].worldPos, neigh.worldPos, Color.yellow);
        }
        */

        if (Input.GetKeyDown(KeyCode.LeftArrow))
        redoAstar(false);
        if (Input.GetKeyDown(KeyCode.RightArrow))
            redoAstar(true);
    }
}
