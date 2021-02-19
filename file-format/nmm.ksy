meta:
  id: nmm
  file-extension: nmm
  endian: le
seq:
  - id: blocks
    type: block
    repeat: until
    repeat-until: _.block_type == block_type::end_of_file
enums:
  state_type:
    0x00: mute
    0x01: normal
    0x02: table
    0x10: codon
    0x11: frame
  block_type:
    0x00: imm_profile
    0x10: nmm_profile
    0xff: end_of_file
types:
  block:
    seq:
      - id: block_type
        type: u1
        enum: block_type
      - id: body
        type:
          switch-on: block_type
          cases:
            'block_type::imm_profile': imm_profile
            'block_type::nmm_profile': nmm_profile
  imm_profile:
    seq:
      - id: abc
        type: abc
      - id: nmodels
        type: u1
      - id: model
        type: model
        repeat: expr
        repeat-expr: nmodels
  nmm_profile:
    seq:
      - id: abc
        type: abc
      - id: nbasep
        type: u2
      - id: basep
        type: basep
        repeat: expr
        repeat-expr: nbasep
      - id: ncodonp
        type: u2
      - id: codonp
        type: codonp
        repeat: expr
        repeat-expr: ncodonp
      - id: ncodonm
        type: u2
      - id: codonm
        type: codonm
        repeat: expr
        repeat-expr: ncodonm
      - id: nmodels
        type: u1
      - id: model
        type: model
        repeat: expr
        repeat-expr: nmodels
  model:
    seq:
      - id: hmm
        type: hmm
      - id: dp
        type: dp
  hmm:
    seq:
      - id: nstates
        type: u2
      - id: states
        type: state
        repeat: expr
        repeat-expr: nstates
      - id: ntransitions
        type: u2
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
        type: f4
      - id: state_type
        type: u1
        enum: state_type
      - id: name_length
        type: u1
      - id: name
        type: str
        encoding: ASCII
        size: name_length + 1
        terminator: 0
      - id: body
        type:
          switch-on: state_type
          cases:
            'state_type::mute': mute_state
            'state_type::normal': normal_state
            'state_type::table': table_state
            'state_type::codon': codon_state
            'state_type::frame': frame_state
  mute_state: {}
  normal_state:
    seq:
      - id: lprobs_size
        type: u1
      - id: lprobs
        type: f4
        repeat: expr
        repeat-expr: lprobs_size
  table_state: {}
  codon_state:
    seq:
      - id: codonp_index
        type: u2
  frame_state:
    seq:
      - id: basep_index
        type: u2
      - id: codonm_index
        type: u2
      - id: epsilon
        type: f4
  transition:
    seq:
      - id: source_state
        type: u2
      - id: target_state
        type: u2
      - id: lprob
        type: f4
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
        type: u2
        repeat: expr
        repeat-expr: max_seq - min_seq + 1
      - id: stride
        type: u2
        repeat: expr
        repeat-expr: max_seq
      - id: size
        type: u2
  dp_emission:
    seq:
      - id: score_size
        type: u4
      - id: score
        type: f4
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
        type: u2
      - id: score
        type: f4
        repeat: expr
        repeat-expr: ntrans
      - id: source_state
        type: u2
        repeat: expr
        repeat-expr: ntrans
      - id: offset_size
        type: u2
      - id: offset
        type: u2
        repeat: expr
        repeat-expr: offset_size
  dp_state_table:
    seq:
      - id: nstates
        type: u2
      - id: min_seq
        type: u1
        repeat: expr
        repeat-expr: nstates
      - id: max_seq
        type: u1
        repeat: expr
        repeat-expr: nstates
      - id: start_lprob
        type: f4
        repeat: expr
        repeat-expr: nstates
      - id: end_state
        type: u2
  array3d:
    seq:
      - id: strides
        type: u2
        repeat: expr
        repeat-expr: 3
      - id: values
        type: f4
        repeat: expr
        repeat-expr: strides[0] * strides[1] * strides[2]
  basep:
    seq:
      - id: lprobs
        type: f4
        repeat: expr
        repeat-expr: 4
  codonp:
    seq:
      - id: lprobs
        type: array3d
  codonm:
    seq:
      - id: symbol_idx
        type: u1
        repeat: expr
        repeat-expr: 128
      - id: lprobs
        type: array3d
