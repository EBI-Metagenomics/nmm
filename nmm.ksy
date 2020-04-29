meta:
  id: nmm
  file-extension: nmm
  endian: le
seq:
  - id: hmm
    type: hmm
  - id: dp
    type: dp
  - id: baset
    type: baset
  - id: codonp
    type: codonp
  - id: codont
    type: codont
enums:
  state_type:
    0: mute
    1: normal
    2: table
    16: codon
    17: frame
types:
  hmm:
    seq:
      - id: abc
        type: abc
      - id: nstates
        type: u4
      - id: states
        type: state
        repeat: expr
        repeat-expr: nstates
      - id: ntransitions
        type: u4
      - id: transitions
        type: transition
        repeat: expr
        repeat-expr: ntransitions
  abc:
    seq:
      - id: abc_type
        type: u1
      - id: nsymbols
        type: u1
      - id: symbols
        type: str
        encoding: ASCII
        size: nsymbols + 1
        terminator: 0
      - id: any_symbol
        type: str
        encoding: ASCII
        size: 1
  state:
    seq:
      - id: start_lprob
        type: f8
      - id: state_type
        type: u1
        enum: state_type
      - id: name_length
        type: u2
      - id: name
        type: str
        encoding: ASCII
        size: name_length + 1
        terminator: 0
      - id: mute_state
        type: mute_state
        if: state_type == state_type::mute
      - id: normal_state
        type: normal_state
        if: state_type == state_type::normal
      - id: table_state
        type: table_state
        if: state_type == state_type::table
      - id: codon_state
        type: codon_state
        if: state_type == state_type::codon
      - id: frame_state
        type: frame_state
        if: state_type == state_type::frame
  mute_state: {}
  normal_state:
    seq:
      - id: lprobs_size
        type: u1
      - id: lprobs
        type: f8
        repeat: expr
        repeat-expr: lprobs_size
  table_state: {}
  codon_state:
    seq:
      - id: codonp_index
        type: u4
  frame_state:
    seq:
      - id: baset_index
        type: u4
      - id: codont_index
        type: u4
      - id: epsilon
        type: f8
  transition:
    seq:
      - id: source_state
        type: u4
      - id: target_state
        type: u4
      - id: lprob
        type: f8
  dp:
    seq:
      - id: seq_code
        type: seq_code
      - id: dp_emission
        type: dp_emission
      - id: dp_trans_table
        type: dp_trans_table
      - id: dp_state_table
        type: dp_state_table
  seq_code:
    seq:
      - id: min_seq
        type: u1
      - id: max_seq
        type: u1
      - id: offset
        type: u4
        repeat: expr
        repeat-expr: max_seq - min_seq + 1
      - id: stride
        type: u4
        repeat: expr
        repeat-expr: max_seq
      - id: size
        type: u4
  dp_emission:
    seq:
      - id: score_size
        type: u4
      - id: score
        type: f8
        repeat: expr
        repeat-expr: score_size
      - id: offset_size
        type: u4
      - id: offset
        type: u4
        repeat: expr
        repeat-expr: offset_size
  dp_trans_table:
    seq:
      - id: ntrans
        type: u4
      - id: score
        type: f8
        repeat: expr
        repeat-expr: ntrans
      - id: source_state
        type: u4
        repeat: expr
        repeat-expr: ntrans
      - id: offset_size
        type: u4
      - id: offset
        type: u4
        repeat: expr
        repeat-expr: offset_size
  dp_state_table:
    seq:
      - id: nstates
        type: u4
      - id: min_seq
        type: u1
        repeat: expr
        repeat-expr: nstates
      - id: max_seq
        type: u1
        repeat: expr
        repeat-expr: nstates
      - id: start_lprob
        type: f8
        repeat: expr
        repeat-expr: nstates
      - id: end_state
        type: u4
  array3d:
    seq:
      - id: strides
        type: u2
        repeat: expr
        repeat-expr: 3
      - id: values
        type: f8
        repeat: expr
        repeat-expr: strides[0] * strides[1] * strides[2]
  baset:
    seq:
      - id: lprobs
        type: f8
        repeat: expr
        repeat-expr: 4
  codonp:
    seq:
      - id: lprobs
        type: array3d
  codont:
    seq:
      - id: symbol_idx
        type: u1
        repeat: expr
        repeat-expr: 128
      - id: lprobs
        type: array3d
