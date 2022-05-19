package edu.duke.cs.banjo.utility;

import java.util.ArrayList;

/**
 * Created by kevindamazyn on 3/25/15.
 */
public class SLUGene {

    int id;
    String name;
    int score;
    ArrayList<SLUGene> parents;

    public SLUGene(int id, String name, int score) {
        this.id = id;
        this.name = name;
        this.score = score;
        this.parents = new ArrayList<SLUGene>();
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public ArrayList<SLUGene> getParents() {
        return parents;
    }

    public void setParents(ArrayList<SLUGene> parents) {
        this.parents = parents;
    }

}
