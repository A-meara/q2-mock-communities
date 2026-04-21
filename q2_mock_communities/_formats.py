import qiime2.plugin.model as model


class CommunityMapFormat(model.TextFileFormat):
    """TSV with sample-id index and a 'community' column (string label)."""

    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
        cols = [c.strip() for c in header.split('\t')]
        if 'community' not in cols:
            raise model.ValidationError(
                "CommunityMapFormat requires a 'community' column."
            )


CommunityMapDirFmt = model.SingleFileDirectoryFormat(
    'CommunityMapDirFmt', 'community_map.tsv', CommunityMapFormat
)
