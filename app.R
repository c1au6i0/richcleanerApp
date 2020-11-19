library(shiny)
library(richcleaner)
library(shinydashboard)
library(shinydashboardPlus)
library(dashboardthemes)
library(shinyFiles)
library(fs)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(dplyr)


jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(
    title = shinyDashboardLogo(
      theme = "blue_gradient",
      boldText = "RichCleaner",
      mainText = "App",
      badgeText = "dev"
    )
  ),
  dashboardSidebar(
    shinyDirButton("directory", "Select Folder", "select_folder"),
    sidebarMenu(
      menuItem("tables", tabName = "tables", icon = icon("table")),
      menuItem("plots",
        tabName = "plots", icon = icon("chart-bar"), startExpanded = TRUE,
        menuSubItem("heatmaps", tabName = "heatmaps"),
        menuSubItem("other graphs", tabName = "rest")
      )
    )
  ),
  dashboardBody(
    # tags$head(tags$style(
    #   HTML("
    #      #buttons {
    #       display: flex;
    #       # align-items: center;
    #       # justify-content: center;
    #       vertical-align: middle;
    #      }
    #      ")
    # )),

    ### changing theme
    shinyDashboardThemes(
      theme = "grey_dark"
    ),
    tabItems(
      tabItem(
        "tables",
        withSpinner(
          DTOutput("tsd_sf")
        )
      ),

      tabItem(
        "heatmaps",
        fluidRow(

          box(
            width = 12,
            column(4, selectInput("value",
              label = "Value to plot",
              choices = c("fdr_q_val", "nes", "nom_p_val", "n_logp_sign"),
              multiple = FALSE
            )),
            column(4, selectInput("gs",
              label = "Gene Set Database",
              choices = "",
              multiple = FALSE
            )),
            column(
              4,
              sliderInput("fdr",
                label = "FDR (log10) threshold",
                min = -10,
                max = -2,
                value = -2
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            column(
              4,
              sliderInput("fontsize_row",
                label = "Fontsize Rows",
                min = 1,
                max = 50,
                value = 12
              )
            ),
            column(
              4,
              sliderInput("fontsize_col",
                label = "Fontsize Columns",
                min = 1,
                max = 50,
                value = 12
              )
            ),
            column(4, style = "margin-top: 25px;",
                   actionButton("plot",
                                  label = "plot",
                                  width = "100%"

            )),
          )
        ),
        fluidRow(
          box(
            width = 12,
            plotOutput("heatmap")

          ))

      ),
      tabItem(
        "rest",
        fluidRow(
          imageOutput("myImage"),
          box(
          HTML("<spam style=color:red;>I am working on this too</spam>"))
        ))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # # get data ----
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = FALSE)

  # reactive values -------------
  dat <- reactiveVal(NULL) # imported data

  observeEvent(
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (length(input$directory) > 1) {
        showModal(modalDialog(
          title = "Alert",
          "Importing data...please wait!",
          footer = NULL,
          size = "m",
          easyClose = FALSE
        ))

        dat_rich <- rich_aggregate(path = parseDirPath(volumes, input$directory))
        dat(dat_rich)
        removeModal()
      }
    }
  )

  output$tsd_sf <- renderDT({
    datatable(dat(),
      selection = "none", filter = "top",
      options = list(autoWidth = TRUE, scrollX = TRUE), width = "auto"
    )
  })


  observe({
    updateSelectInput(session, "gs", choices = unique(dat()$gs))
  })

  pheatmap <- eventReactive(input$plot,{

    rich_df <- as.data.frame(req(dat()))

    rich_plot  <- tryCatch(
          rich_pheatmap(dat = rich_df,
                        fdr_threshold = 10^input$fdr,
                        gs = input$gs,
                        value = input$value,
                        fontsize_row = input$fontsize_row,
                        fontsize_col= input$fontsize_col),
                       error = function(c) conditionMessage(c)
      )

    if (rich_plot == "must have n >= 2 objects to cluster"){
      to_plot <-  ggplot() +
                theme_void() +
                geom_text(aes(0 , 0, label = "No enough clusters to plot with these filtering parameters!"), size = 10) +
                xlab(NULL)

    } else {
      to_plot <-   rich_plot
    }

    to_plot
  })


  output$heatmap <- renderPlot({
    pheatmap()
  })

  output$myImage <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- normalizePath("/Users/heverz/Desktop/GSEA_PreRanked_EpitheliumOnly/Enterocytes_dkonormal/Biocarta/Enterocytes_dkonormal_Biocarta.GseaPreranked.1605217181475/enplot_BIOCARTA_MTOR_PATHWAY_99.png")

    # Return a list containing the filename and alt text
    list(src = filename,
         alt = "x")

  }, deleteFile = FALSE)



}

# Run the application
shinyApp(ui = ui, server = server)
