library(ggrepel)
if (F) {
    p + geom_text_repel(aes(x, y, label, max.overlaps = 20))
    p + geom_label_repel(aes(x, y, label, max.overlaps = 20))
}
