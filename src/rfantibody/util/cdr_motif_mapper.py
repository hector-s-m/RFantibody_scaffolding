"""
CDR motif mapper for aromatic peptide motif scaffolding.

Maps motif residues to positions within a CDR loop, computes flanking
residue counts from the sampled CDR loop length, and generates
fixed-position dicts for ProteinMPNN.
"""

import random

import numpy as np


class CDRMotifMapper:
    """
    Maps a peptide motif into a CDR loop within an AbPose.

    The flanking residue counts are derived from the CDR loop length
    (already set by adjust_loop_lengths) minus the motif length. The
    motif is centered in the loop with a random ±1 offset.
    """

    def __init__(self, pose, cdr_loop: str, motif_data: dict):
        """
        Args:
            pose: AbPose object (after adjust_loop_lengths has been called)
            cdr_loop: CDR loop name, e.g. 'H3'
            motif_data: dict from motif_loader.load_motif_combined_pdb()
        """
        self.pose = pose
        self.cdr_loop = cdr_loop.upper()
        self.motif_data = motif_data
        self.motif_len = len(motif_data['motif_seq'])

        # Will be set by place_motif_in_loop()
        self.flank_n = None
        self.flank_c = None
        self.motif_global_indices = None

    def place_motif_in_loop(self) -> tuple:
        """
        Compute flank sizes from the CDR loop length and motif length.

        The CDR loop length has already been sampled by adjust_loop_lengths().
        Flanks are split roughly evenly with a random ±1 offset.

        Returns:
            (flank_n, flank_c): number of designed residues before and after motif

        Raises:
            ValueError: if the CDR loop is too short for the motif
        """
        loop_start, loop_end = self.get_cdr_loop_positions()
        loop_len = loop_end - loop_start

        flank_total = loop_len - self.motif_len
        if flank_total < 4:
            raise ValueError(
                f"CDR loop {self.cdr_loop} has {loop_len} residues, too short "
                f"for motif of {self.motif_len} residues (need at least "
                f"{self.motif_len + 4} for minimum 2 flanking residues on each side). "
                f"Increase the loop length range in design_loops."
            )

        # Split evenly with random offset
        half = flank_total // 2
        remainder = flank_total % 2
        if remainder > 0:
            # Randomly assign the extra residue to N or C flank
            offset = random.choice([0, 1])
        else:
            offset = 0

        self.flank_n = half + offset
        self.flank_c = flank_total - self.flank_n

        # Compute global indices of the motif within the full pose
        self.motif_global_indices = list(
            range(loop_start + self.flank_n,
                  loop_start + self.flank_n + self.motif_len)
        )

        print(f"Motif scaffolding: {self.cdr_loop} loop length={loop_len}, "
              f"motif length={self.motif_len}, "
              f"flanks=({self.flank_n}, {self.flank_c})")

        return self.flank_n, self.flank_c

    def get_cdr_loop_positions(self) -> tuple:
        """
        Get the (start, end) global indices of the CDR loop in the full pose.

        Returns:
            (start, end): start is inclusive, end is exclusive
        """
        # Get the combined loop masks from the pose
        loop_masks = self.pose.combine_loop_masks()

        if self.cdr_loop not in loop_masks:
            raise ValueError(
                f"Unknown CDR loop: {self.cdr_loop}. "
                f"Valid options: {list(loop_masks.keys())}"
            )

        mask = loop_masks[self.cdr_loop]
        indices = np.where(mask)[0]

        if len(indices) == 0:
            raise ValueError(
                f"CDR loop {self.cdr_loop} has no residues in the pose"
            )

        return int(indices[0]), int(indices[-1]) + 1

    def get_fixed_positions_for_mpnn(self) -> dict:
        """
        Generate fixed position dict for ProteinMPNN.

        Returns positions of motif residues within their chain, using
        1-indexed residue numbers relative to chain start (matching
        the convention in SampleFeatures.loop_string2fixed_res).

        Returns:
            dict: {chain_letter: [1-indexed positions to fix]}
        """
        if self.motif_global_indices is None:
            raise RuntimeError("Call place_motif_in_loop() first")

        # Determine which chain the CDR belongs to
        chain = self.cdr_loop[0]  # 'H' or 'L'

        # The motif global indices are relative to the full pose (H+L+T).
        # For ProteinMPNN, we need chain-relative 1-indexed positions.
        if chain == 'H':
            chain_offset = 0
        elif chain == 'L':
            chain_offset = self.pose.H.seq.shape[0] if self.pose.has_H() else 0
        else:
            raise ValueError(f"Motif CDR must be on H or L chain, got {chain}")

        fixed_positions = [
            (idx - chain_offset) + 1  # Convert to 1-indexed chain-relative
            for idx in self.motif_global_indices
        ]

        return {chain: fixed_positions}
