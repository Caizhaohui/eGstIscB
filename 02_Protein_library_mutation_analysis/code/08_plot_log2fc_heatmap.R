#' ==============================================================================
#' Protein library 深度突变扫描 (DMS) 热图无缝优化可视化脚本 (左侧生化分类与紧凑全景版)
#' 
#' 更新与优化要点:
#' 1. Y轴氨基酸分类与标记移至热图左侧 (参考图示规范):
#'    - Small: A, G
#'    - Hydrophob.: M, V, L, I
#'    - Nuc.: T, S, C
#'    - Amide: Q, N
#'    - Arom.: Y, W, F
#'    - Neg.: E, D
#'    - Pos.: R, K, H
#'    - 未分类/特殊: P, * (终止密码子，数据编码为 X)
#'    - 左侧配置分类名 (竖排 90度) + 归类分组折角括号 + 氨基酸单字母标签 (紧贴热图左侧)
#'    - 热图右侧完全移除冗余行名标签，保持视觉极致清爽与对齐
#' 2. 单元格高度设置为 cellheight = 10 (cellwidth = 5.5)
#' 3. 横坐标以 50 间隔标记刻度:
#'    - Part 1 (Pos 2 ~ 213): 2, 50, 100, 150, 200
#'    - Part 2 (Pos 214 ~ 424): 从 214 开始以 50 为步长: 214, 264, 314, 364, 414
#' 4. 双行堆叠全景图 (Combined) 视觉紧凑优化:
#'    - 上下两行间距紧贴 (gap_pt = 10，彻底消除空白断层)
#'    - 全景图不显示主标题，仅保留单一颜色标签条 (右上侧)
#'    - 全景图全部文字使用 Arial，并以 cairo_pdf 嵌入 TrueType，便于 Illustrator 编辑
#'    - 左右边界与行名垂直方向绝对对齐
#' ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(grid)
  library(gtable)
})

#' 读取 FASTA 参考序列
read_reference_sequence <- function(fasta_path = "reference.fasta") {
  if (!file.exists(fasta_path)) {
    stop(paste("未找到参考序列文件:", fasta_path))
  }
  fasta_lines <- readLines(fasta_path, warn = FALSE)
  seq_str <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
  return(seq_str)
}

#' 获取指定位点的野生型氨基酸
get_wt_amino_acid <- function(pos, ref_seq) {
  idx <- pos + as.integer(Sys.getenv("REFERENCE_OFFSET", "6"))  # Historical default; set REFERENCE_OFFSET for another mapping.
  if (idx >= 1 && idx <= nchar(ref_seq)) {
    return(substr(ref_seq, idx, idx))
  }
  return(NA_character_)
}

#' 构建符合生化分类排序与连续位点的热图矩阵
build_heatmap_matrix <- function(data_clean, ref_seq, pos_start, pos_end, aa_order) {
  grid_df <- expand_grid(
    Position = pos_start:pos_end,
    Mut_AA = aa_order
  ) %>%
    rowwise() %>%
    mutate(WT_AA = get_wt_amino_acid(Position, ref_seq)) %>%
    ungroup() %>%
    mutate(is_wildtype = (WT_AA == Mut_AA))
  
  mat_df <- grid_df %>%
    left_join(
      data_clean %>% select(Position, Mut_AA, log2_fc_plot),
      by = c("Position", "Mut_AA")
    ) %>%
    mutate(log2_fc_plot = ifelse(is_wildtype, NA_real_, log2_fc_plot)) %>%
    select(Position, Mut_AA, log2_fc_plot) %>%
    pivot_wider(names_from = Position, values_from = log2_fc_plot) %>%
    column_to_rownames("Mut_AA")
  
  mat_df <- mat_df[aa_order, ]
  return(as.matrix(mat_df))
}

#' 创建热图左侧生化分类与氨基酸标识 Grob
#' @param show_aa 是否绘制 AA 单字母标签
#' @param show_cat 是否绘制分类名文字；FALSE 时只保留折角括号
#' @param show_bracket 是否绘制折角括号（竖线 + tick）；FALSE 时整体留白
create_left_annotation_grob <- function(fontsize_aa = 6.5, fontsize_cat = 6.5,
                                        fontfamily = NULL, show_aa = TRUE,
                                        show_cat = TRUE, show_bracket = TRUE) {
  text_gp <- function(fontsize) {
    gp <- gpar(fontsize = fontsize, fontface = "plain", col = "black")
    if (!is.null(fontfamily)) {
      gp$fontfamily <- fontfamily
    }
    gp
  }
  # 21 种突变类型展示标签 (自上而下)
  aa_display <- c(
    "A", "G",             # Small
    "M", "V", "L", "I",   # Hydrophob.
    "T", "S", "C",        # Nuc.
    "Q", "N",             # Amide
    "Y", "W", "F",        # Arom.
    "E", "D",             # Neg.
    "R", "K", "H",        # Pos.
    "P",                  # 特殊
    "*"                   # 终止密码子
  )
  n_rows <- length(aa_display)
  
  # 分类区间定义 (1-based index)
  # 无 AA 字母且 cell 较矮时用略缩写名，避免竖排文字挤出分组带
  if (isTRUE(show_aa)) {
    categories <- list(
      list(name = "Small",      start = 1,  end = 2),
      list(name = "Hydrophob.", start = 3,  end = 6),
      list(name = "Nuc.",       start = 7,  end = 9),
      list(name = "Amide",      start = 10, end = 11),
      list(name = "Arom.",      start = 12, end = 14),
      list(name = "Neg.",       start = 15, end = 16),
      list(name = "Pos.",       start = 17, end = 19)
    )
  } else {
    categories <- list(
      list(name = "Small",   start = 1,  end = 2),
      list(name = "Hydroph.", start = 3,  end = 6),
      list(name = "Nuc.",    start = 7,  end = 9),
      list(name = "Amide",   start = 10, end = 11),
      list(name = "Arom.",   start = 12, end = 14),
      list(name = "Neg.",    start = 15, end = 16),
      list(name = "Pos.",    start = 17, end = 19)
    )
  }
  
  # 计算在 0 到 1 归一化视口中的 y 坐标 (0 为底，1 为顶)
  get_cell_center <- function(i) (n_rows - i + 0.5) / n_rows
  get_cell_top    <- function(i) (n_rows - i + 1.0) / n_rows
  get_cell_bottom <- function(i) (n_rows - i) / n_rows
  
  grobs <- list()
  
  # 1. 氨基酸单字母标签 (可选；紧贴热图左边缘)
  if (isTRUE(show_aa)) {
    for (i in seq_along(aa_display)) {
      grobs[[paste0("aa_", i)]] <- textGrob(
        label = aa_display[i],
        x = unit(0.85, "npc"),
        y = unit(get_cell_center(i), "npc"),
        just = c("center", "centre"),
        gp = text_gp(fontsize_aa)
      )
    }
  }
  
  # 组间微小间距 (避免连续两组线条粘连)
  # cellheight=5.5 时相对原 10pt 等比缩小；无 AA 字母时略收紧
  gap_npc <- if (isTRUE(show_aa)) 1.2 / 210 else 0.7 / 210
  
  # 布局坐标：有 AA 字母时括号靠左；无 AA 时整体右移贴紧热图
  if (isTRUE(show_aa)) {
    x_vline <- 0.58
    x_tick_end <- 0.72
    x_text <- 0.30
  } else {
    x_vline <- 0.72
    x_tick_end <- 0.92
    x_text <- 0.38
  }
  
  # 2. 分类折角括号与分类文本 (竖排 90度)
  for (cat in categories) {
    y_top <- get_cell_top(cat$start) - gap_npc
    y_bot <- get_cell_bottom(cat$end) + gap_npc
    y_center <- (y_top + y_bot) / 2
    
    # 竖向连线 + 顶部横向折角 tick（可整体关闭）
    if (isTRUE(show_bracket)) {
      grobs[[paste0("line_", cat$name)]] <- linesGrob(
        x = unit(c(x_vline, x_vline), "npc"),
        y = unit(c(y_bot, y_top), "npc"),
        gp = gpar(col = "black", lwd = 0.8)
      )
      grobs[[paste0("tick_", cat$name)]] <- linesGrob(
        x = unit(c(x_vline, x_tick_end), "npc"),
        y = unit(c(y_top, y_top), "npc"),
        gp = gpar(col = "black", lwd = 0.8)
      )
    }
    
    # 分类名称文本 (竖向逆时针旋转 90 度，自下而上阅读)
    if (isTRUE(show_cat)) {
      grobs[[paste0("text_", cat$name)]] <- textGrob(
        label = cat$name,
        x = unit(x_text, "npc"),
        y = unit(y_center, "npc"),
        rot = 90,
        just = c("center", "centre"),
        gp = text_gp(fontsize_cat)
      )
    }
  }
  
  return(gTree(
    children = do.call(gList, grobs),
    name = "left_annotation",
    vp = viewport(clip = "off")
  ))
}

#' 将左侧分类 Grob 挂载到 pheatmap 的 gtable 中
attach_left_annotation <- function(p, width_pt = 48, fontfamily = NULL,
                                   show_aa = TRUE, fontsize_aa = 6.5,
                                   fontsize_cat = 6.5, show_cat = TRUE,
                                   show_bracket = TRUE) {
  gt <- p$gtable
  # 设定 column 2 宽度为 width_pt (约 0.67 吋 / 无 AA 时约 0.42 吋)
  gt$widths[2] <- unit(width_pt, "bigpts")
  anno_grob <- create_left_annotation_grob(
    fontsize_aa = fontsize_aa,
    fontsize_cat = fontsize_cat,
    fontfamily = fontfamily,
    show_aa = show_aa,
    show_cat = show_cat,
    show_bracket = show_bracket
  )
  gt <- gtable_add_grob(
    gt,
    grobs = anno_grob,
    t = 4, l = 2, b = 4, r = 2,
    z = Inf,
    name = "left_annotation",
    clip = "off"
  )
  return(gt)
}

#' 生成横坐标刻度标签
create_col_labels <- function(positions, mode = c("part1", "part2_seq50", "part2_round50")) {
  mode <- match.arg(mode)
  if (mode == "part1") {
    # 2, 50, 100, 150, 200
    labels <- ifelse(positions %% 50 == 0 | positions == min(positions), as.character(positions), "")
  } else if (mode == "part2_seq50") {
    # 从 214 开始，每 50 间隔: 214, 264, 314, 364, 414
    labels <- ifelse((positions - min(positions)) %% 50 == 0, as.character(positions), "")
  } else if (mode == "part2_round50") {
    # 214 + 50 的倍数: 214, 250, 300, 350, 400
    labels <- ifelse(positions %% 50 == 0 | positions == min(positions), as.character(positions), "")
  }
  return(labels)
}

#' 将横坐标位点数字改为竖排，并在对应列中心加刻度
#' 相当于水平文字顺时针旋转 90 度：数字开头靠近刻度，自上而下阅读
add_vertical_col_axis <- function(gt, labels, fontsize = 10, fontfamily = "Arial",
                                  tick_pt = 9, gap_pt = 3, pad_pt = 4) {
  labels <- as.character(labels)
  idx <- which(nzchar(labels))
  if (length(idx) == 0) {
    return(gt)
  }
  n <- length(labels)
  x <- (idx - 0.5) / n
  labs <- labels[idx]
  longest <- labs[which.max(nchar(labs))]
  span_pt <- suppressWarnings(convertWidth(
    grobWidth(textGrob(
      longest,
      gp = gpar(fontsize = fontsize, fontfamily = fontfamily)
    )),
    "bigpts",
    valueOnly = TRUE
  ))

  ticks <- segmentsGrob(
    x0 = unit(x, "npc"),
    x1 = unit(x, "npc"),
    y0 = unit(1, "npc"),
    y1 = unit(1, "npc") - unit(tick_pt, "bigpts"),
    gp = gpar(col = "black", lwd = 1.15, lineend = "butt")
  )
  txt <- textGrob(
    label = labs,
    x = unit(x, "npc"),
    y = unit(1, "npc") - unit(tick_pt + gap_pt, "bigpts"),
    just = c("right", "centre"),
    rot = 90,  # 水平文字顺时针转 90 度，自上而下阅读
    gp = gpar(
      fontsize = fontsize,
      fontfamily = fontfamily,
      fontface = "plain",
      col = "black"
    )
  )
  axis <- gTree(
    children = gList(ticks, txt),
    name = "vertical_col_axis",
    vp = viewport(
      layout.pos.row = 5,
      layout.pos.col = 3,
      clip = "off",
      height = unit(tick_pt + gap_pt + span_pt + pad_pt, "bigpts"),
      just = "bottom"
    )
  )

  # show_colnames = FALSE 时 col_names 行被压到 0pt，把刻度层挂在矩阵行向下伸出，
  # 这样热图高度和页面高度都保持不变。
  hit <- which(gt$layout$name == "matrix")
  if (length(hit) == 1) {
    gt$grobs[[hit]] <- gTree(
      children = gList(gt$grobs[[hit]], axis),
      name = "matrix_with_axis"
    )
  } else {
    hit <- which(gt$layout$name == "col_names")
    if (length(hit) != 1) {
      stop("未找到矩阵或 col_names，无法添加竖排刻度")
    }
    gt$grobs[[hit]] <- axis
  }
  gt
}

#' 覆写 pheatmap 矩阵 cell 边框线宽与颜色
#' pheatmap 不支持直接设 cell 线宽，需修改其内部 rect grob 的 gp
#' 注意：rect 的 gp$lwd 默认为 NULL（走 grid 默认 0.75），直接 $lwd <- 会被丢弃，
#' 必须显式构建 gpar 覆写。
set_matrix_border_lwd <- function(gt, lwd = 0.35, col = "grey60") {
  hit <- which(gt$layout$name %in% c("matrix", "matrix_with_axis", "matrix_with_separators"))
  if (length(hit) != 1) {
    stop("未找到矩阵，无法调整边框线宽")
  }
  m <- gt$grobs[[hit]]
  # 向下钻取找到 rect 子对象并就地替换
  replace_rect_gp <- function(g) {
    if (inherits(g, "rect")) {
      g$gp <- do.call(
        gpar,
        utils::modifyList(unclass(g$gp), list(lwd = lwd, col = col))
      )
      class(g$gp) <- "gpar"
      return(g)
    }
    if (!is.null(g$children)) {
      g$children <- do.call(gList, lapply(g$children, replace_rect_gp))
    }
    g
  }
  gt$grobs[[hit]] <- replace_rect_gp(m)
  gt
}

#' 在热图矩阵区叠加生化分类间水平分隔线
#' 边界位于第 2|3、6|7、9|10、11|12、14|15、16|17 行之间
add_category_separator_lines <- function(gt, n_rows = 21, lwd = 0.6, col = "black") {
  # 分类结束行 (1-based, 自上而下)；P 与 * 不分组，不画线
  after_rows <- c(2, 6, 9, 11, 14, 16)
  # grid y: 0 为底；第 i 行底边 = (n_rows - i) / n_rows
  y_npc <- (n_rows - after_rows) / n_rows
  seps <- segmentsGrob(
    x0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y0 = unit(y_npc, "npc"),
    y1 = unit(y_npc, "npc"),
    gp = gpar(col = col, lwd = lwd, lineend = "butt")
  )
  sep_tree <- gTree(children = gList(seps), name = "category_separators")

  hit <- which(gt$layout$name %in% c("matrix", "matrix_with_axis"))
  if (length(hit) != 1) {
    stop("未找到矩阵，无法添加分类分隔线")
  }
  gt$grobs[[hit]] <- gTree(
    children = gList(gt$grobs[[hit]], sep_tree),
    name = "matrix_with_separators"
  )
  gt
}

#' 水平位点轴：数字水平排布 + 短刻度线（CCW90 版式）
#' 锚定在矩阵底边向下伸出 (clip=off)，避免叠在热图中间
add_horizontal_col_axis <- function(gt, labels, fontsize = 16, fontfamily = "Arial",
                                    tick_pt = 9, gap_pt = 3, pad_pt = 4) {
  labels <- as.character(labels)
  idx <- which(nzchar(labels))
  if (length(idx) == 0) {
    return(gt)
  }
  n <- length(labels)
  x <- (idx - 0.5) / n
  labs <- labels[idx]
  span_pt <- suppressWarnings(convertHeight(
    grobHeight(textGrob(
      "8",
      gp = gpar(fontsize = fontsize, fontfamily = fontfamily)
    )),
    "bigpts",
    valueOnly = TRUE
  ))
  axis_h <- tick_pt + gap_pt + span_pt + pad_pt

  ticks <- segmentsGrob(
    x0 = unit(x, "npc"),
    x1 = unit(x, "npc"),
    y0 = unit(1, "npc"),
    y1 = unit(1, "npc") - unit(tick_pt, "bigpts"),
    gp = gpar(col = "black", lwd = 1.15, lineend = "butt")
  )
  txt <- textGrob(
    label = labs,
    x = unit(x, "npc"),
    y = unit(1, "npc") - unit(tick_pt + gap_pt, "bigpts"),
    just = c("center", "top"),
    rot = 0,
    gp = gpar(
      fontsize = fontsize,
      fontfamily = fontfamily,
      fontface = "plain",
      col = "black"
    )
  )
  # y=0 + just=top: 视口顶边贴矩阵底边，内容向下画出矩阵外
  axis <- gTree(
    children = gList(ticks, txt),
    name = "horizontal_col_axis",
    vp = viewport(
      x = unit(0.5, "npc"),
      y = unit(0, "npc"),
      width = unit(1, "npc"),
      height = unit(axis_h, "bigpts"),
      just = c("centre", "top"),
      clip = "off"
    )
  )

  hit <- which(gt$layout$name %in% c("matrix", "matrix_with_axis", "matrix_with_separators"))
  if (length(hit) == 1) {
    gt$grobs[[hit]] <- gTree(
      children = gList(gt$grobs[[hit]], axis),
      name = "matrix_with_axis",
      vp = viewport(clip = "off")
    )
    # gtable 默认裁切矩阵格，必须关闭才能画出底边外的刻度
    if ("clip" %in% names(gt$layout)) {
      gt$layout$clip[hit] <- "off"
    }
  } else {
    hit <- which(gt$layout$name == "col_names")
    if (length(hit) != 1) {
      stop("未找到矩阵或 col_names，无法添加水平刻度")
    }
    gt$heights[gt$layout$t[hit]] <- unit(axis_h, "bigpts")
    gt$grobs[[hit]] <- axis
  }
  gt
}

#' 将顶层 gtable 中的图例整体下移 offset_in 英寸（水平位置不变）
shift_legend_down <- function(gt, offset_in) {
  hit <- which(gt$layout$name == "legend")
  if (length(hit) != 1 || offset_in <= 0) {
    return(gt)
  }
  leg <- gt$grobs[[hit]]
  gt$grobs[[hit]] <- gTree(
    children = gList(leg),
    vp = viewport(
      y = unit(0.5, "npc") - unit(offset_in, "inches"),
      just = "centre",
      clip = "off"
    ),
    name = "legend_shifted"
  )
  if ("clip" %in% names(gt$layout)) {
    gt$layout$clip[hit] <- "off"
  }
  gt
}

#' 精确绘制紧凑堆叠双行全景热图
#' @param legend_shift_down_in 顶层图例下移量（英寸）；使图例落到两行热图竖直中线
draw_stacked_panoramic_heatmap <- function(gt_top, gt_bottom, gap_pt = 10,
                                           legend_shift_down_in = 0) {
  grid.newpage()
  
  h_top_in <- convertHeight(sum(gt_top$heights), "in", valueOnly = TRUE)
  h_bot_in <- convertHeight(sum(gt_bottom$heights), "in", valueOnly = TRUE)
  gap_in   <- gap_pt / 72
  total_h_in <- h_top_in + gap_in + h_bot_in

  # 下移量：顶层中心 → 整体中心，等于 (h_bot + gap) / 2
  if (legend_shift_down_in > 0) {
    gt_top <- shift_legend_down(gt_top, legend_shift_down_in)
  }
  
  # 全局居中容器
  vp_container <- viewport(
    x = 0.5, y = 0.5,
    width = unit(1, "npc"),
    height = unit(total_h_in, "inches")
  )
  pushViewport(vp_container)
  
  # 顶层热图；clip=off 允许图例画到面板下方
  vp_top <- viewport(
    x = 0.5,
    y = unit(h_bot_in + gap_in + h_top_in / 2, "inches"),
    width = unit(1, "npc"),
    height = unit(h_top_in, "inches"),
    clip = "off"
  )
  pushViewport(vp_top)
  grid.draw(gt_top)
  popViewport()
  
  # 底层热图
  vp_bot <- viewport(
    x = 0.5,
    y = unit(h_bot_in / 2, "inches"),
    width = unit(1, "npc"),
    height = unit(h_bot_in, "inches")
  )
  pushViewport(vp_bot)
  grid.draw(gt_bottom)
  popViewport()
  
  popViewport()
  return(invisible(total_h_in))
}

#' 执行优化主流程
run_pheatmap_pipeline <- function(
    input_csv = "results/log2_foldchange.csv",
    fasta_file = "reference.fasta",
    out_dir = "figures_optimized",
    condition_name = "comparison",
    part2_tick_mode = "part2_seq50",
    gap_pt = 10) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  message(">>> 正在读取参考序列: ", fasta_file)
  ref_seq <- read_reference_sequence(fasta_file)
  
  message(">>> 正在读取与清洗突变数据: ", input_csv)
  data_raw <- read_csv(input_csv, show_col_types = FALSE)
  
  data_clean <- data_raw %>%
    mutate(
      WT_AA = str_sub(Mutation, 1, 1),
      Mut_AA = str_sub(Mutation, -1, -1),
      Position = as.numeric(str_extract(Mutation, "(?<=^[A-Z])-?\\d+(?=[A-Z\\*]$)"))
    ) %>%
    filter(!is.na(Position), Position >= 2, Position <= 424) %>%
    mutate(
      is_wildtype = WT_AA == Mut_AA,
      log2_fc_plot = ifelse(is_wildtype, NA_real_, log2_foldchange)
    )
  
  # 按照图片分类严格定义的氨基酸顺序 (自上而下)
  aa_data_order <- c(
    "A", "G",             # Small
    "M", "V", "L", "I",   # Hydrophob.
    "T", "S", "C",        # Nuc.
    "Q", "N",             # Amide
    "Y", "W", "F",        # Arom.
    "E", "D",             # Neg.
    "R", "K", "H",        # Pos.
    "P",                  # 特殊
    "X"                   # 终止密码子 (数据编码为 X)
  )
  
  message(">>> 正在构建 Part 1 (Pos 2~213) 与 Part 2 (Pos 214~424) 数据矩阵...")
  mat1 <- build_heatmap_matrix(data_clean, ref_seq, 2, 213, aa_data_order)
  mat2 <- build_heatmap_matrix(data_clean, ref_seq, 214, 424, aa_data_order)
  
  # 统一全局色阶与 breaks (严格沿用 filterBylog0.5 规范)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_fc <- max(c(as.vector(mat1), as.vector(mat2)), na.rm = TRUE)
  
  breaks <- c(
    seq(-3, 0.5, length.out = 50),
    seq(0.5 + 1e-6, max_fc, length.out = 50)
  )
  breaks <- unique(breaks)
  message(sprintf(">>> 色阶范围: [-3.0, %.3f], 临界点: 0.5, breaks 阶数: %d", max_fc, length(breaks)))
  
  # 横坐标 50 个间隔刻度标签
  pos_vec1 <- as.numeric(colnames(mat1))
  labels1 <- create_col_labels(pos_vec1, mode = "part1")
  
  pos_vec2 <- as.numeric(colnames(mat2))
  labels2 <- create_col_labels(pos_vec2, mode = part2_tick_mode)
  
  message(">>> Part 1 横坐标刻度标记: ", paste(labels1[labels1 != ""], collapse = ", "))
  message(">>> Part 2 横坐标刻度标记: ", paste(labels2[labels2 != ""], collapse = ", "))
  
  # ----------------------------------------------------------------------------
  # 1. 生成独立 Part 1 热图 (Pos 2 ~ 213, cellheight = 10, 左侧生化分类，独立图例)
  # ----------------------------------------------------------------------------
  message(">>> 正在渲染独立 Part 1 (Positions 2 to 213)...")
  p1_indep <- pheatmap(
    mat1,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = paste0(condition_name, " - Part 1 (Positions 2 to 213)"),
    treeheight_row = 0,
    cellwidth = 5.5,
    cellheight = 10,
    fontsize_col = 6,
    fontsize = 9,
    labels_col = labels1,
    breaks = breaks,
    show_rownames = FALSE,  # 隐藏右侧行名
    legend = TRUE,
    silent = TRUE
  )
  gt1_indep <- attach_left_annotation(p1_indep, width_pt = 48)
  
  part1_pdf <- file.path(out_dir, paste0(condition_name, "_Part1_Pos2-213_filterBylog0.5.pdf"))
  part1_png <- file.path(out_dir, paste0(condition_name, "_Part1_Pos2-213_filterBylog0.5.png"))
  
  pdf(part1_pdf, width = 20, height = 4.5)
  grid::grid.draw(gt1_indep)
  invisible(dev.off())
  
  png(part1_png, width = 20, height = 4.5, units = "in", res = 300)
  grid::grid.draw(gt1_indep)
  invisible(dev.off())
  
  message("  [已保存] ", part1_pdf)
  message("  [已保存] ", part1_png)
  
  # ----------------------------------------------------------------------------
  # 2. 生成独立 Part 2 热图 (Pos 214 ~ 424, cellheight = 10, 左侧生化分类，独立图例)
  # ----------------------------------------------------------------------------
  message(">>> 正在渲染独立 Part 2 (Positions 214 to 424)...")
  p2_indep <- pheatmap(
    mat2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = paste0(condition_name, " - Part 2 (Positions 214 to 424)"),
    treeheight_row = 0,
    cellwidth = 5.5,
    cellheight = 10,
    fontsize_col = 6,
    fontsize = 9,
    labels_col = labels2,
    breaks = breaks,
    show_rownames = FALSE,  # 隐藏右侧行名
    legend = TRUE,
    silent = TRUE
  )
  gt2_indep <- attach_left_annotation(p2_indep, width_pt = 48)
  
  part2_pdf <- file.path(out_dir, paste0(condition_name, "_Part2_Pos214-424_filterBylog0.5.pdf"))
  part2_png <- file.path(out_dir, paste0(condition_name, "_Part2_Pos214-424_filterBylog0.5.png"))
  
  pdf(part2_pdf, width = 20, height = 4.5)
  grid::grid.draw(gt2_indep)
  invisible(dev.off())
  
  png(part2_png, width = 20, height = 4.5, units = "in", res = 300)
  grid::grid.draw(gt2_indep)
  invisible(dev.off())
  
  message("  [已保存] ", part2_pdf)
  message("  [已保存] ", part2_png)
  
  # ----------------------------------------------------------------------------
  # 3. 生成双行堆叠全景拼图 (Combined):
  #    - 顶层 p1_comb: 无主标题，显示唯一颜色标签条，挂载左侧分类
  #    - 底层 p2_comb: 无标题，不显示颜色标签条，挂载左侧分类
  #    - 全部文字为 Arial（cairo_pdf 嵌入 TrueType，Illustrator 中可直接编辑）
  #    - 两图左右列宽完全对齐，纵向紧凑贴合 (gap_pt = 10)
  # ----------------------------------------------------------------------------
  message(">>> 正在渲染双行紧凑全景图专用组件 (无标题、Arial、唯一图例)...")
  
  # pheatmap 在空 pdf 设备上测量文字，该设备没有 Arial，会回退到 Helvetica 度量并给出警告。
  # 文字对象本身仍记录为 Arial，最终由 cairo_pdf 按 Arial 绘制并嵌入。
  p1_comb <- suppressWarnings(pheatmap(
    mat1,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = NA,                       # 全景图不显示标题
    treeheight_row = 0,
    cellwidth = 5.5,
    cellheight = 10,
    fontsize_col = 6,
    fontsize = 10,
    labels_col = labels1,
    breaks = breaks,
    show_rownames = FALSE,           # 隐藏右侧行名
    show_colnames = FALSE,           # 位点数字由竖排刻度统一绘制
    legend = TRUE,                   # 顶层显示唯一颜色标签条
    fontfamily = "Arial",
    silent = TRUE
  ))
  
  p2_comb <- suppressWarnings(pheatmap(
    mat2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = NA,                       # 底层无标题
    treeheight_row = 0,
    cellwidth = 5.5,
    cellheight = 10,
    fontsize_col = 6,
    fontsize = 10,
    labels_col = labels2,
    breaks = breaks,
    show_rownames = FALSE,           # 隐藏右侧行名
    show_colnames = FALSE,           # 位点数字由竖排刻度统一绘制
    legend = FALSE,                  # 底层关闭图例
    fontfamily = "Arial",
    silent = TRUE
  ))
  
  gt_top <- attach_left_annotation(p1_comb, width_pt = 48, fontfamily = "Arial")
  gt_bot <- attach_left_annotation(p2_comb, width_pt = 48, fontfamily = "Arial")
  
  # 强制使底层 p2 的第 5 列宽度与顶层 p1 的图例列宽度一致，实现左右绝对对齐
  gt_bot$widths[5] <- gt_top$widths[5]
  
  comb_pdf <- file.path(out_dir, paste0(condition_name, "_Combined_Pos2-424_filterBylog0.5.pdf"))
  comb_png <- file.path(out_dir, paste0(condition_name, "_Combined_Pos2-424_filterBylog0.5.png"))
  
  # 精确计算两图加上间距的物理高度
  # 空 pdf 设备无法识别 Arial，度量时会警告；此处压掉该警告，绘制仍走 cairo。
  h_top_in <- suppressWarnings(convertHeight(sum(gt_top$heights), "in", valueOnly = TRUE))
  h_bot_in <- suppressWarnings(convertHeight(sum(gt_bot$heights), "in", valueOnly = TRUE))
  gap_in   <- gap_pt / 72
  canvas_h_in <- round(h_top_in + gap_in + h_bot_in + 0.3, 2)
  message(sprintf(">>> 双行堆叠全景图画布高度自动设定为: %.2f 吋 (top: %.2f, gap: %.2f, bot: %.2f)",
                  canvas_h_in, h_top_in, gap_in, h_bot_in))
  
  # cairo_pdf 将 Arial 作为 TrueType 文字嵌入，Illustrator 中可直接点选编辑
  suppressWarnings({
    cairo_pdf(comb_pdf, width = 20, height = canvas_h_in, family = "Arial")
    draw_stacked_panoramic_heatmap(gt_top, gt_bot, gap_pt = gap_pt)
    invisible(dev.off())

    png(comb_png, width = 20, height = canvas_h_in, units = "in", res = 300,
        type = "cairo", family = "Arial")
    draw_stacked_panoramic_heatmap(gt_top, gt_bot, gap_pt = gap_pt)
    invisible(dev.off())
  })
  
  message("  [已保存] ", comb_pdf)
  message("  [已保存] ", comb_png)

  # --------------------------------------------------------------------------
  # 3b. 竖排刻度版本：不覆盖上面的全景图
  #     位点数字字号 10，顺时针旋转 90 度，并在该列加刻度
  # --------------------------------------------------------------------------
  message(">>> 正在渲染竖排刻度全景图 (新文件，不覆盖原图)...")
  axis_fontsize <- 10
  gt_top_v <- add_vertical_col_axis(
    gt_top, labels1,
    fontsize = axis_fontsize, fontfamily = "Arial"
  )
  gt_bot_v <- add_vertical_col_axis(
    gt_bot, labels2,
    fontsize = axis_fontsize, fontfamily = "Arial"
  )

  h_top_v <- suppressWarnings(convertHeight(sum(gt_top_v$heights), "in", valueOnly = TRUE))
  h_bot_v <- suppressWarnings(convertHeight(sum(gt_bot_v$heights), "in", valueOnly = TRUE))
  canvas_v_in <- round(h_top_v + gap_in + h_bot_v + 0.3, 2)
  message(sprintf(">>> 竖排刻度全景图画布高度: %.2f 吋", canvas_v_in))

  comb_pdf_v <- file.path(
    out_dir,
    paste0(condition_name, "_Combined_Pos2-424_filterBylog0.5_verticalAxis.pdf")
  )
  suppressWarnings({
    cairo_pdf(comb_pdf_v, width = 20, height = canvas_v_in, family = "Arial")
    draw_stacked_panoramic_heatmap(gt_top_v, gt_bot_v, gap_pt = gap_pt)
    invisible(dev.off())
  })
  message("  [已保存] ", comb_pdf_v)

  message(">>> 全部热图渲染更新成功！")
  return(invisible(list(
    gt1_indep = gt1_indep,
    gt2_indep = gt2_indep,
    gt_top = gt_top,
    gt_bot = gt_bot,
    canvas_height = canvas_h_in
  )))
}

#' CCW90 优化版全景热图管线
#' - cell 5.5 x 5.5 pt 正方形
#' - 去每格描边，仅在生化分类间加 6 条水平分隔线
#' - 水平位点轴 16 pt，左侧删除 AA 单字母（仅保留分类名）
#' - 输出新文件，不覆盖原图
run_ccw90_optimized_pipeline <- function(
    input_csv = "results/log2_foldchange.csv",
    fasta_file = "reference.fasta",
    out_dir = "figures_optimized",
    condition_name = "comparison",
    part2_tick_mode = "part2_seq50",
    gap_pt = 10,
    cell_size = 5.5,
    axis_fontsize = 16,
    border_lwd = 0.1,
    left_width_pt = 4) {

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  message(">>> [CCW90优化] 读取参考序列: ", fasta_file)
  ref_seq <- read_reference_sequence(fasta_file)

  message(">>> [CCW90优化] 读取与清洗突变数据: ", input_csv)
  data_raw <- read_csv(input_csv, show_col_types = FALSE)

  data_clean <- data_raw %>%
    mutate(
      WT_AA = str_sub(Mutation, 1, 1),
      Mut_AA = str_sub(Mutation, -1, -1),
      Position = as.numeric(str_extract(Mutation, "(?<=^[A-Z])-?\\d+(?=[A-Z\\*]$)"))
    ) %>%
    filter(!is.na(Position), Position >= 2, Position <= 424) %>%
    mutate(
      is_wildtype = WT_AA == Mut_AA,
      log2_fc_plot = ifelse(is_wildtype, NA_real_, log2_foldchange)
    )

  aa_data_order <- c(
    "A", "G",
    "M", "V", "L", "I",
    "T", "S", "C",
    "Q", "N",
    "Y", "W", "F",
    "E", "D",
    "R", "K", "H",
    "P",
    "X"
  )

  message(">>> [CCW90优化] 构建矩阵 Part1(2~213) / Part2(214~424)...")
  mat1 <- build_heatmap_matrix(data_clean, ref_seq, 2, 213, aa_data_order)
  mat2 <- build_heatmap_matrix(data_clean, ref_seq, 214, 424, aa_data_order)

  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_fc <- max(c(as.vector(mat1), as.vector(mat2)), na.rm = TRUE)
  breaks <- unique(c(
    seq(-3, 0.5, length.out = 50),
    seq(0.5 + 1e-6, max_fc, length.out = 50)
  ))
  message(sprintf(">>> 色阶范围: [-3.0, %.3f], 临界点: 0.5", max_fc))

  pos_vec1 <- as.numeric(colnames(mat1))
  labels1 <- create_col_labels(pos_vec1, mode = "part1")
  pos_vec2 <- as.numeric(colnames(mat2))
  labels2 <- create_col_labels(pos_vec2, mode = part2_tick_mode)
  message(">>> Part 1 刻度: ", paste(labels1[labels1 != ""], collapse = ", "))
  message(">>> Part 2 刻度: ", paste(labels2[labels2 != ""], collapse = ", "))

  message(">>> [CCW90优化] 渲染 pheatmap (cell=", cell_size,
          "pt, border_lwd=", border_lwd, ", 无左侧文字, 位点轴=", axis_fontsize, "pt)...")
  p1_comb <- suppressWarnings(pheatmap(
    mat1,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = NA,
    treeheight_row = 0,
    cellwidth = cell_size,
    cellheight = cell_size,
    fontsize_col = 6,
    fontsize = 10,
    labels_col = labels1,
    breaks = breaks,
    show_rownames = FALSE,
    show_colnames = FALSE,
    legend = TRUE,
    fontfamily = "Arial",
    silent = TRUE
  ))

  p2_comb <- suppressWarnings(pheatmap(
    mat2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    na_col = "white",
    border_color = "grey60",
    main = NA,
    treeheight_row = 0,
    cellwidth = cell_size,
    cellheight = cell_size,
    fontsize_col = 6,
    fontsize = 10,
    labels_col = labels2,
    breaks = breaks,
    show_rownames = FALSE,
    show_colnames = FALSE,
    legend = FALSE,
    fontfamily = "Arial",
    silent = TRUE
  ))

  # 删除全部左侧标签与括号（AA 字母 + 分类名 + 折角），左列收窄
  gt_top <- attach_left_annotation(
    p1_comb,
    width_pt = left_width_pt,
    fontfamily = "Arial",
    show_aa = FALSE,
    show_cat = FALSE,
    show_bracket = FALSE
  )
  gt_bot <- attach_left_annotation(
    p2_comb,
    width_pt = left_width_pt,
    fontfamily = "Arial",
    show_aa = FALSE,
    show_cat = FALSE,
    show_bracket = FALSE
  )
  gt_bot$widths[5] <- gt_top$widths[5]

  # cell 边框改细（覆写矩阵 rect 的 gp$lwd）
  gt_top <- set_matrix_border_lwd(gt_top, lwd = border_lwd, col = "grey60")
  gt_bot <- set_matrix_border_lwd(gt_bot, lwd = border_lwd, col = "grey60")

  gt_top <- add_horizontal_col_axis(
    gt_top, labels1,
    fontsize = axis_fontsize, fontfamily = "Arial"
  )
  gt_bot <- add_horizontal_col_axis(
    gt_bot, labels2,
    fontsize = axis_fontsize, fontfamily = "Arial"
  )

  # 两行热图净间隙目标 28 pt。
  # 实测：净间隙 ≈ stack_gap_pt + 10 pt（gtable 在矩阵底边之下还有约 10 pt）。
  # 故 stack_gap_pt = 18 时，热图底边到下一行顶边约为 28 pt，刚好容纳 16 pt 位点数字。
  axis_reserve_pt <- 18
  stack_gap_pt <- axis_reserve_pt
  gap_in <- stack_gap_pt / 72
  h_top <- suppressWarnings(convertHeight(sum(gt_top$heights), "in", valueOnly = TRUE))
  h_bot <- suppressWarnings(convertHeight(sum(gt_bot$heights), "in", valueOnly = TRUE))
  # 图例下移到两行整体竖直中线：偏移 = (底层高 + 间距) / 2
  legend_shift_down_in <- (h_bot + gap_in) / 2
  # 底层轴伸出画布底部；16pt 数字会再向下约 10pt，单独加底边留白，不改变两行间距
  canvas_h_in <- round(h_top + gap_in + h_bot + axis_reserve_pt / 72 + 0.25 + 22 / 72, 2)
  message(sprintf(
    ">>> 画布高度: %.2f 吋 (top: %.2f, gap: %.2f, bot: %.2f, legend_down: %.2f)",
    canvas_h_in, h_top, gap_in, h_bot, legend_shift_down_in
  ))

  out_pdf <- file.path(
    out_dir,
    paste0(condition_name, "_Combined_Pos2-424_filterBylog0.5_verticalAxis_CCW90_optimized.pdf")
  )
  out_png <- file.path(
    out_dir,
    paste0(condition_name, "_Combined_Pos2-424_filterBylog0.5_verticalAxis_CCW90_optimized.png")
  )

  suppressWarnings({
    cairo_pdf(out_pdf, width = 20, height = canvas_h_in, family = "Arial")
    draw_stacked_panoramic_heatmap(
      gt_top, gt_bot,
      gap_pt = stack_gap_pt,
      legend_shift_down_in = legend_shift_down_in
    )
    invisible(dev.off())

    png(out_png, width = 20, height = canvas_h_in, units = "in", res = 300,
        type = "cairo", family = "Arial")
    draw_stacked_panoramic_heatmap(
      gt_top, gt_bot,
      gap_pt = stack_gap_pt,
      legend_shift_down_in = legend_shift_down_in
    )
    invisible(dev.off())
  })

  message("  [已保存] ", out_pdf)
  message("  [已保存] ", out_png)
  message(">>> [CCW90优化] 完成！")
  return(invisible(list(
    gt_top = gt_top,
    gt_bot = gt_bot,
    canvas_height = canvas_h_in,
    pdf = out_pdf,
    png = out_png
  )))
}

# ------------------------------------------------------------------------------
# 主执行入口
# ------------------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3 || length(args) > 4) {
    stop("Usage: Rscript 08_plot_log2fc_heatmap.R INPUT_CSV REFERENCE_FASTA OUTPUT_DIR [COMPARISON_LABEL]")
  }
  label <- if (length(args) == 4) args[[4]] else "comparison"
  run_ccw90_optimized_pipeline(input_csv = args[[1]], fasta_file = args[[2]],
                               out_dir = args[[3]], condition_name = label,
                               part2_tick_mode = "part2_seq50", gap_pt = 10)
}
