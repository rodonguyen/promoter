package qut;

public class TaskHandler {
    private final Gene referenceGene;
    private final Gene gene;
    private final GenbankRecord record;

    public TaskHandler(Gene referenceGene, Gene gene, GenbankRecord record) {
        this.referenceGene = referenceGene;
        this.gene = gene;
        this.record = record;
    }

    public Gene getReferenceGene() {
        return referenceGene;
    }
    public Gene getGene() {
        return gene;
    }
    public GenbankRecord getRecord() {
        return record;
    }
}
